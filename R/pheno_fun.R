# Phenotyping functions ---------------

#' This script contains functions with phenotyping utilities. That is:
#'  - Reading ancestral alleles by linking parental and NAM populations
#'  - Creating effects based on a series of models
#'  - Applying these effects to obtain heretability-controlled phenotypes and rescaled
#'    effects.
#'
#'  Developed by Alejandro Thérèse Navarro
#'  October 2019

link_NAM <- function(#function to connect ancestral founder alleles and NAM alleles
  crossfile, #founderallele file of the cross
  marker=NULL, #index of marker or markers to get
  totallele, #path to Total_pop.txt, that is a pedigreeSim-type table with chromosomes and individuals, specifying alleles for each parental chromosome
  parental=F, #F:return the whole cross NAM, T:returns only the parental alleles
  parents=10, #number of parents
  ploidy=4
){
  #Obtain chromosome names of parents
  headcross <- unlist(data.table::fread(crossfile,nrow=1,header=F)) #first column is "marker"
  par_names <- headcross[(1+(1:(parents*ploidy)))]#first columns must contain parent names
  markers <- data.table::fread(crossfile,select="marker",header=T)

  #find genotypes in datafile with true alleles
  par_founder <- data.table::fread(totallele,select=par_names,header=T)
  par_founder <- as.matrix(par_founder)

  if(!parental){
    #read the cross genotypes
    cross_geno <- data.table::fread(crossfile,header=T)

    cross_geno <- as.matrix(cross_geno[,-1,with=F]) #turn into matrix for speed
    #use chromosome numbers in cross_geno as indeces to on par_founder, which contains the true alleles ordered so that they coincide with parental chromosomes
    result <- t(sapply(1:nrow(markers),function(x) par_founder[x,cross_geno[x,]+1]))
    rownames(result) <- unlist(markers)
    colnames(result) <- headcross[-1] #change names of matrix, exclude the "marker column"
  }else if(parental){
    rownames(par_founder) <- marker
    return(par_founder)
  }else{stop("Parentals must be either T or F")}

  if(!is.null(marker)) result <- result[marker,] #if certain markers are selected, only return those

  return(result)
}

#' Effect generator
#'
#' For each group of alleles belonging to an ancestral group, we can define a method
#' to sample effects. Currently two methods are supported, "anc_average" in which
#' alleles from each ancestral group are sampled from a normal distribution around
#' an ancestral average (with spread equal to anc_sd); and "fix_eff" in which
#' each ancestral group has `n` non-zero alleles, of effect `size` (each ancestral
#' group's allele size is specified as a vector).
#'
#' @param gen genotypes per chromosome per individual
#' @param anc_alleles matrix (per column)/list defining which
#' alleles belong to which ancestral group
#' @param anc_sd numeric, indicating sd for the anc_average method
#' @param seed
#' @param method "anc_average" or "fix_eff"
#' @param size num or numeric vector, identifying the size of the active allele
#' per ancestral group. It is recycled.
#' @param n number of active alleles per ancestral group.
#'
#' @return
#' @export
#'
#' @examples
#' gen <- matrix(c(1,1,2,3,5,4,3,2,2,3,4,1,3,3,1),nrow=1)
#' multi.effects(gen,
#'               anc_alleles = matrix(1:12,ncol=4), method="fix_eff",
#'               size = c(3,-2),n=1)
effect_gen <- function(
  gen, #genotypes (per chromosome of individual)
  anc_alleles, #matrix (per column)/list defining which alleles belong to which ancestral group
  anc_sd=1,
  seed=7,
  method = "anc_average", #or "fix_eff"
  size = 3, #size of the effect per ancestral group
  n = 1 #number of alleles with an effect
){

  if(is.matrix(anc_alleles)) anc_alleles <- mat2list(anc_alleles)
  if(is.vector(gen)) gen <- matrix(gen,nrow=1)

  effects <- lapply(1:nrow(gen),function(k){
    #unique alleles present in locus k
    a <- unique(gen[k,])

    #to which ancestral does each allele correspond
    anc_matrix <- sapply(anc_alleles,function(c) a%in%c)
    rownames(anc_matrix) <- a

    #Number of alleles present in which ancestral groups
    anc_n <- colSums(anc_matrix)[colSums(anc_matrix)!=0]

    # Average per anc group routine -----
    #Define the ancestral averages
    if(method == "anc_average"){
      set.seed(seed*k)
      anc_mean <- runif(length(anc_n),0,1)

      eff <- lapply(1:length(anc_n),function(i){
        set.seed(seed*k*i)
        #randomly generate effects with a specific seed
        e <- rnorm(anc_n[i],anc_mean[i],sd=anc_sd)
        #which group are we taking (to put the right allele names)
        group <- which(colSums(anc_matrix)!=0)[i]
        names(e) <- rownames(anc_matrix)[which(anc_matrix[,group])]
        return(e)
      })
      names(eff) <- names(anc_n)

    }else if(method == "fix_eff"){
      #Specific size per allele ------
      if(length(anc_n) > length(size)) size <- rep(size,length(anc_n))
      eff <- lapply(1:length(anc_n), function(i){
        e <- c(rep(size[i],n),rep(0,anc_n[i]-n))

        #which group are we taking (to put the right allele names)
        group <- which(colSums(anc_matrix)!=0)[i]
        names(e) <- rownames(anc_matrix)[which(anc_matrix[,group])]
        return(e)
      })
      names(eff) <- names(anc_n)
    }


    return(do.call(c,eff))
  })

  return(effects)
}


pheno<-function(
  genotypes, #alleles per chromosome (loci x chromosome matrix)
  effects, #as many effects as rows has the genotype matrix
  ploidy,
  #The polygen effect is no longer considered useful
  polygen=NULL, #additional polygenic genotypes
  polygen.effects=NULL, #polygenic effects can be specified
  polygen.sd=mu/50/nrow(polygen), #variation of polygenic effects
  seed=7,
  herit=0.5,
  mu=100,
  Evar=2,
  return.effects=F
){
  if(is.vector(genotypes)) genotypes<-matrix(genotypes,nrow=1,
                                             dimnames = list(NULL,names(genotypes)))
  if(!is.list(effects)) effects<-list(effects)

  #Phenotype without environmental variance
  phenotype<-phsum(genotypes,effects,ploidy)

  #After calculating phenotype, we would like to adapt it
  #so that it has the heritability and the mean we want
  set.seed(seed)
  env<-rnorm(length(phenotype),mean=0,sd=Evar)
  Sg<-var(phenotype)
  Se<-var(env)

  #Calculation of heritability control parameter
  b2<-Se*herit/(Sg*(1-herit))#calculation of b squared
  effects2<-lapply(effects,function(x) x*sqrt(b2))#rescaling of genetic effects
  phenotype<-phsum(genotypes,effects2,ploidy)

  #Calculation of the mean control parameter
  c<-(mu-mean(phenotype))/(ploidy*nrow(genotypes))
  effects2<-lapply(effects2,function(x) x+c)#rescaling of genetic effects
  phenotype<-phsum(genotypes,effects2,ploidy)

  if(!is.null(polygen)){

    #Define polygenic effects if they are not defined
    if(is.null(polygen.effects)){
      effects.polyg<-apply(polygen,1,function(gen){
        #mu should be zero, sd affects how much different alleles are
        set.seed(seed)
        e<-rnorm(length(unique(gen)),0,polygen.sd)
        names(e)<-unique(gen)
        return(e)
      })
    }

    #apply polygenic effects to calculate polygenic term
    polygenic<-phsum(polygen,polygen.effects,ploidy)
  }

  #add up genetic and environmental effects
  result<-phenotype+env
  if(!is.null(polygen)) result<-result+polygenic

  if(return.effects){
    return(list(pheno=result,
                effects=effects2))
  }
  return(result)
}

phsum <- function(
  genotypes,
  effects,
  ploidy,
  table=F #whether to give the results in table format (one for each locus) or already summed
){
  if(is.vector(genotypes)) genotypes<-matrix(genotypes,nrow=1,
                                             dimnames = list(NULL,names(genotypes)))
  if(!is.list(effects)) effects<-list(effects)

  phenolist<-lapply(1:nrow(genotypes),function(k){
    gen<-genotypes[k,]
    ef<-effects[[k]]

    #incidence matrix (allele per chromosome)
    am<-sapply(gen,function(g) g==unique(gen))
    rownames(am)<-unique(gen)

    #reorder to match the order of the effects (must have same names)
    am<-t(am)[ ,names(ef)]
    cef<-am%*%ef #chromosome effect vector

    #Calculate phenotype at this locus (sum chromosomes per individual)
    ph<-sapply(1:(length(gen)/ploidy), function(i){
      j<-1:ploidy+(i-1)*ploidy
      r<-sum(cef[j]) #additive model (sums effects at each chromosome)

      names(r)<-paste(Reduce(
        #find the common characters in the names of each group of chromosomes
        intersect2, strsplit(names(cef[j,]), NULL)),
        collapse = '')

      return(r)
    })
    return(ph)
  })
  if(table){
    return(Reduce(cbind,phenolist))
  }else{
    return(Reduce(`+`,phenolist))
  }

}

# Misc ---------------------------------

#' These functions are not phenotyping functions but they help
#' the functions above

mat2list <- function(matrix,rowise=F){
  if(rowise){
    #if rowise we transpose the matrix so that rows are columns
    matrix<-t(matrix)
  }

  n <- ncol(matrix)
  res <- lapply(1:n,function(i) matrix[,i])

  return(res)
}

intersect2<-function (x, y)
{
  y <- as.vector(y)
  y[match(as.vector(x), y, 0L)]
}
