# PedigreeSim handling functions ----------------

#' PedigreeSim is a software for generating polyploid meiosis simulations
#' and it requires a series of files to be run. In order to call the
#' program using R-style programming, a series of functions have
#' been developed.
#'
#' Developed by Alejandro Therese Navarro
#' October 2019

makedir <- function(folder){if(!dir.exists(folder)){dir.create(folder)}}

#' Make PedSim Parameter File
#'
#' @param output adress of the output to be created (without extension),
#' @param chrom adress of the chromosome file. If not specified, it will assume the same name as output.
#' @param pedfile adress of the pedigree (.ped) file. If not specified, it will assume the same name as output.
#' @param mapgen adress of the .map file. If not specified, it will assume the same name as output.
#' @param founder adress of the .gen file. If not specified, it will assume the same name as output.
#' @param outname adress of the output files. If not specified, it will assume the same name as output.
#' @param ploidy ploidy level to be used
#' @param mapfun mapping function to be used. Either "HALDANE" or "KOSAMBI"
#' @param miss character indivating missing data
#' @param pairing pairing type to be used. 1 if probabilistic, 0 if determined by .chrom file.
#'
#' @return creates a parameter for PedigreeSim.
#' @export
#'
#' @examples
makePar <- function(output,chrom=output,pedfile=output,mapgen=output,founder=output,
                    outname=output,ploidy=4,mapfun="HALDANE", miss="NA",pairing=1){
  #Writes the file
  write(file=paste0(output,".par"),
        c(paste0("PLOIDY = ",ploidy),
          paste0("MAPFUNCTION = ",mapfun),
          paste0("MISSING = ",miss),
          paste0("CHROMFILE = ",chrom,".chrom"),
          paste0("PEDFILE = ",pedfile,".ped"),
          paste0("MAPFILE = ",mapgen,".map"),
          paste0("FOUNDERFILE = ",founder,".gen"),
          paste0("OUTPUT = ",outname),
          paste0("NATURALPAIRING = ",pairing)))
}

# Pedigree file generation --------------

#' Make NAM pedigree file
#' @description Generates a NAM pedigree file. Names of the offspring are coded as
#' O(ffspring)number_crossnumber.
#'
#' @param parents character vector of parent names
#' @param offspring number of offpsring per NAM
#' @param output character string for pedigree adress
#' @param prefix prefix to add to offspring names
#'
#' @return A .ped file of NAM structure for PedigreeSim
#' @export
#'
#' @examples
makeNAMPed<- function(
  parents,
  offspring,
  output="Pedigree",
  prefix=NULL
){

  #Create the rows regarding parents
  pedPar <- data.frame(name = parents, parent1 = NA, parent2 = NA)

  if(length(offspring) < length(parents[-1])) offspring <- rep(offspring,length(parents))

  #Create the rows regarding the offspring. Each offspring will be a cross of P1 with the rest
  pedNAM <- lapply(2:length(parents),function(i) {
    names <- paste0(prefix,"C",i-1,"_","O",sprintf("%02.f",1:offspring[i-1]))
    data.frame(name = names,parent1 = parents[1],parent2 = parents[i],
               stringsAsFactors = F)
  })
  pedNAM <- do.call(rbind,pedNAM)
  pedNAM <- rbind(pedPar,pedNAM)

  write.table(pedNAM,paste0(output,".ped"),quote = F,row.names = F, col.names = T)

}

#' Make Random Pedigree File
#' @description Generates a pedigree file of a random mixing population. A set of founders
#' are created labelled as G0# (generation 0), and they are randomly paired, without selfing,
#' to produce the following generations.
#'
#' @param output name of the file to be created
#' @param nfounder number of founder individuals
#' @param popsize population size of each generation
#' @param ngener nubmer of generations to be made
#'
#' @return A PedigreeSim pedigree file following a random mating process.
#' @export
#'
#' @examples
makeRanPed<- function(
  output,
  nfounder=10,
  popsize=100,
  ngener=10
){
  #lists to store the different generations
  Name <- list(); Parent1 <- list(); Parent2 <- list()

  #Create names for as many founders as specified (P1, P2...)
  Name[[1]] <- paste0("G0_",sprintf("%02d",1:nfounder))
  #create empty lists with slots for as many founders as specified
  Parent1[[1]] <- rep(NA, nfounder)
  Parent2[[1]] <- Parent1[[1]]

  #loop that generates offsprings for each generation, in Name the individual names are stored.
  #Name[1] is a vector with parents, Name[2] is a vector with first generation (G1), etc.
  #Parent1 and 2 are sampled from previous generation, with replacement
  for (g in 2:(ngener+1)) {
    Name[[g]] <- paste0("G", g-1, "_", sprintf("%03d",1:popsize))
    Parent1[[g]] <- sample(Name[[g-1]], popsize, replace=TRUE)
    Parent2[[g]] <- sample(Name[[g-1]], popsize, replace=TRUE)
    #to prevent selfings, if Parent1 and parent2 are the same, parent 2 is resampled.
    repeat {
      self <- Parent1[[g]] == Parent2[[g]] #self=vector with T or F if they are the same.
      if (sum(self) == 0) break() #if all are F, sum(self)==0
      Parent2[[g]][self] <- sample(Name[[g-1]], sum(self), replace=TRUE) #only those elements with self=true are resampled
    }
  }
  #a dataframe is made that contains the pedigree
  ped <- data.frame(Name=unlist(Name), Parent1=unlist(Parent1), Parent2=unlist(Parent2))
  #the df is written into the output file
  write.table(ped,paste0(output,".ped"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
}

makeAncPed<-function(
  output,
  nfounder=10,
  popsize=100,
  ngener=10,
  groups=1,
  geneflow=0,
  seed=7
){
  #ask giorgio's opinion
  #We divide by two because there are two parents
  flowprob <- (groups-1)*geneflow/groups/2

  #Naming for the ancestral groups
  if(groups>length(LETTERS)){
    agname<-unlist(
      sapply(LETTERS,function(x) paste0(x,LETTERS))
    )
  }else{
    agname<-LETTERS
  }

  #List that will store, for each generation, names of individuals
  Name<-list(lapply(1:groups,function(i){
    ag<-agname[i]
    paste0(ag,"0_",sprintf("%02d",1:nfounder))
  }))
  founders<-unlist(Name)

  #List that stores, for each generation, names of parents (first gen empty)
  Parent1<-list(lapply(1:groups,function(i) rep(NA,nfounder)))
  Parent2<-Parent1

  #loop that generates offsprings for each generation, in Name the individual names are stored.
  #Name[1] is a vector with parents, Name[2] is a vector with first generation (G1), etc.
  #Parent1 and 2 are sampled from previous generation, with replacement
  for (g in 2:(ngener+1)) {
    Name[[g]] <- lapply(1:groups,function(i){
      ag<-agname[i]
      paste0(ag,g-1,"_",sprintf("%03d",1:popsize))
    })

    #GENE FLOW DETERMINATION
    #each individual has a local parent, or a flowed parent (from another AG)
    set.seed(seed*g)
    flow1<-sapply(1:groups,function(x) runif(popsize)<flowprob)

    #to avoid individuals with two flowed parents
    reps<-0
    repeat{
      reps<-reps+1
      set.seed(seed*g*2+reps)
      flow2<-sapply(1:groups,function(x) runif(popsize)<flowprob)
      out<-flow1&flow2
      if(sum(out)==0) break()
    }

    Parent1[[g]]<-lapply(1:groups,function(ag){
      thisag<-Name[[g-1]][[ag]]
      otherag<-unlist(Name[[g-1]][-ag])
      flo<-flow1[,ag]

      set.seed(seed*g)
      pa<-sample(thisag,popsize,replace=TRUE)
      set.seed(seed*g+1)
      pa[flo]<-sample(otherag,sum(flo),replace=T)
      return(pa)
    })

    Parent2[[g]]<-lapply(1:groups,function(ag){
      thisag<-Name[[g-1]][[ag]]
      otherag<-unlist(Name[[g-1]][-ag])
      flo<-flow2[,ag]

      set.seed(seed*g*2)
      pa<-sample(thisag,popsize,replace=T)
      set.seed(seed*g*2+1)
      pa[flo]<-sample(otherag,sum(flo),replace=T)

      #We add this check to prevent selfing
      reps<-0
      repeat{
        reps<-reps+1
        self<-pa==Parent1[[g]][[ag]]
        if(sum(self)==0) break()
        set.seed(seed*g*2+reps)
        pa[self&!flo]<-sample(thisag,sum(self&!flo),replace=T)
        set.seed(seed*g*2+reps^2)
        pa[self&flo]<-sample(otherag,sum(self&flo),replace=T)
      }
      return(pa)
    })

  }
  #a dataframe is made that contains the pedigree
  ped <- data.frame(Name=unlist(Name), Parent1=unlist(Parent1), Parent2=unlist(Parent2))
  #the df is written into the output file
  write.table(ped,paste0(output,".ped"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

  return(founders)
}

# Other file generation ------------------

#' Make a PedigreeSim chromosome file
#'
#' @param output name of the file to be created
#' @param chms vector of chromosome names
#' @param cM numeric vector of chromosome lengths, in centimorgans
#' @param centros vector of centromere positions
#' @param pref preferential pairing value (all chromosomes will have the same value)
#' @param quad quadrivalent fraction (all chromosomes will have the same value)
#'
#' @return A PedigreeSim chromosome file.
#' @export
#'
#' @examples
makeChrom <- function(
  output, #name or adress of the file to be created
  chms, #vector of chromosome names
  cM, #vector of chromosome length
  centros, #vector of centromere positions
  pref, #preferential pairing (all chr the same value)
  quad #quadrivalent fraction (all chromosomes the same)
){

  #check that all vectors have the same length
  if(!length(chms)&&length(cM)&&length(centros)==(length(chms)))
  {stop("Unequal input vector lengths")}

  ## Create header:
  write(file=paste0(output,".chrom"),
        c("chromosome length centromere prefPairing quadrivalents"))

  #print each row with each vector
  if(length(quad) < length(chms)) quad <- rep(quad,length(chms))
  if(length(pref) < length(chms)) pref <- rep(pref,length(chms))

  for(i in 1:length(chms)){
    write(file=file.path(paste0(output,".chrom")),
          paste(chms[i],cM[i],centros[i],pref[i],quad[i]),append=TRUE)
  }
}

#' Make PedigreeSim map file
#'
#' @param output name of the file to be created
#' @param marker vector with marker names
#' @param chr vector with chromosome of each marker
#' @param pos vector with marker positions (in each chromosome)
#' @param head logical value indicating whether the header "marker chromosome position"
#' should be written or not. Useful for writing maps by chromosome, instead of all at once.
#' @param body logical value indicating whether the body (each marker) should be written or not.
#' Useful for writing maps by chromosome, instead of all at once.
#'
#' @return A PedigreeSim genetic map file
#' @export
#'
#' @examples
makeMap <- function(
  output, #name of file to be written
  marker, #vector with marker names
  chr, #vector with the chr of each marker
  pos, #vector with marker positions
  head=T, #write the head ("marker chromosome position") or not
  body=T #logical value for whether to write the body or not
){

  #Create header. depends on head.
  if(head){write(file=paste0(output,".map"),
                 c("marker chromosome position"))}

  #writes rows with marker name, chr to which they belong and position. Depends on body.
  if(body){
    for(i in 1:length(marker)){
      write(file=paste0(output,".map"),
            paste(marker[i],chr[i],pos[i]),append=TRUE)
    }
  }
}

#' Make a PedigreeSim .gen file
#' @description Formats a data.frame into a .gen file as required by PedigreeSim. The data.frame passed
#' needs to have a first column named "markers" which contains marker names.
#'
#' @param output name of the file to be created
#' @param data data frame where each row is a marker (first column must be named "markers" and contain
#' marker names) and each column is a parental homologue, #' labelled as Parent_1,
#' Parent_2 ... Parent_ploidy. Each column will contain an allele state (numbers,
#' letters, whatever code).
#'
#' @return A PedigreeSim gen file
#' @export
#'
#' @examples
makeGen<- function(
  output, #output file name
  data #data frame. Each row is a marker. First column is named "markers" and has marker names. Rest of columns contain the alleles of each parental homologue, labelled such as Parent_1, Parent_2 etc..
){

  #check right object file cause otherwise write.table will not work
  if(!is.data.frame(data)){print("data file is not dataframe"); break}

  #needed for PedigreeSim, checks that first column is called markers (important header in the .gen file)
  names<-colnames(data)
  if(names[1]!="marker"){print("First column of data must be named 'markers', check column names"); break}

  write.table(file=paste0(output,".gen"),data,quote=F,row.names=F)
}

# Random parameters -------------

#' Random PedigreeSim chromosome file
#' @description For testing purposes. Will generate a random chromosome file with
#' specific chromosome number, random lengths of chromosomes, optional random centromere positions,
#' with specified preferential and quadrivalent probabilities.
#'
#' @param output name of the file to be created
#' @param chms number of chromosomes to be simulated
#' @param CM average chromosome length
#' @param sdCM standard deviation to be applied to chromosome length
#' @param centr average or fixed centromere position (2 is 1/2, 3 is 1/3 and so on)
#' @param rancentr logical value indicating whether centromere position should be
#' random, or fixed.
#' @param pref probability of preferential pairing.
#' @param quad probability of quadrivalent pairing.
#'
#' @return A PedigreeSim chromosome file with random chromosome structure.
#' @export
#'
#' @examples
ranChrom <- function(
  output, #name of the file
  chms=1, #number of chromosomes
  CM=100, #average chr length
  sdCM=CM/2, #st dev of chr length
  centr=2, #average centromere position (2=half, 3= third and so on) or centromere position (in cM) if rancentr=F
  rancentr=T, #(T/F) whether centr are random
  pref=0, #probability of preferential pairing,
  quad=0 #probability of quadrivalents
){


  #each chromosome has a letter name. If there are more than 26 chromosomes, double letter names will be used (AA, AB...)
  lets<-c()
  for(i in 1:length(LETTERS)){
    lets<-c(lets,paste0(LETTERS[i],LETTERS))
  }
  chms<-seq(1:chms)
  if(length(chms)>26){chms<-lets[chms]}else{ chms<-LETTERS[chms]}


  #make vector with chr length. If random, with mean CM and sd=sdCM, otherwise just all same length
  repeat{
    chrlength<-c()
    chrlength<-c(rnorm(length(chms),mean=CM,sd=sdCM))
    if(all(chrlength>0)){break}
  }
  #order CM according to length
  chrlength<-chrlength[order(chrlength,decreasing = T)]

  #Make random centromeres vector. Centros are always on the top half of the chromosome length,
  #and are prevented to be negative. If random, mean=chr length/centr sd=chr length/6
  centros<-c()
  if(rancentr){
    for(i in 1:length(chms)){
      repeat{
        centros[i]<-rnorm(1,mean=chrlength[i]/centr,sd=chrlength[i]/6)
        if(centros[i]>0){break}
      }
      if((centros[i]-chrlength[i]/2)>0){centros[i]<-chrlength[i]-centros[i]}
    }
  }else{centros<-rep(centr,length(chms))}

  #make the file
  makeChrom(output=output,chms=chms,cM=chrlength,centros=centros,pref=pref,quad=quad)
}

#' Random PedigreeSim map file
#' @description For testing purposes. A random map file is generated, based on a chromosome file.
#' Markers are randomly position but at a specified density along the chromosomes
#'
#' @param chrom .chrom file adress of the file to base the genetic map on.
#' @param output name of the output file. Defaults to the same name as the chromosome file,
#' with different termination.
#' @param density marker density to simulate, in markers/cM
#'
#' @return A random genetic map file, based on a chromosome file.
#' @export
#'
#' @examples
ranMap<-function(
  chrom, #chrom file name
  output=chrom, #output name
  density=1 #marker density in markers/cM
){
  #read chrom file to get the chromosome names and the chromosome lengths
  table<-read.delim(chrom,sep=" ")

  chromname<-as.vector(table[,1])
  chromlength<-as.vector(table[,2])

  #a file is made, but only the header is written
  makeMap(output,body=F)

  #for each chromosome
  for(i in 1:length(chromname)){
    #we generate an interval vector containing interval division points.
    interval<-seq(from=0,to=chromlength[i],by=1/density)
    marker<-c();markchrom<-c();markname<-c()

    #for each interval, a single marker is positioned randomly (runif), a vector with the chromosome name is made
    #and a marker name is generated such as A001, first marker of chromosome A.
    for(j in 1:(length(interval)-1)){
      marker[j]<-runif(1,interval[j],interval[j+1])
      markchrom[j]<-chromname[i]
      markname[j]<-paste0(chromname[i],sprintf("%03d",j))

    }

    #for each chromosome, the markers are written in the marker file (without header), and the loop continues.
    makeMap(output,marker=markname,chr=markchrom,pos=marker,head=F)
  }
}

#' Make a random PedigreeSim .gen file
#' @description For testing purposes. Based on a genetic map file, a number of parental
#' genotypes is simulated such that the genotypes originate from a flat distribution. Ploidy
#' of the parents can be specified.
#'
#' @param map genetic map to base the simulation on (with extension)
#' @param parents number of parents to be simulated
#' @param output name of the output file. Defaults to same name as map, with different extension
#' @param ploidy ploidy number of the parents
#'
#' @return A random PedigreeSim .gen file for a specified number of parents.
#' @export
#'
#' @examples
ranGen<- function(
  map, #name of the map file to be used (with extension)
  parents, #number of parents
  output=map, #name of the output
  ploidy=4,
  object_map
){
  #read map file to obtain marker names.
  if(object_map){
    markers <- map
  }else{
    markers<-read.delim(map,sep=" ")
  }


  #create data frame where the genotypes will be stored. First column contains marker names.
  gen<-data.frame("marker"=markers[,1])
  #the rest of columns are created, P01_1, P01_2, etc.
  parname<-c()
  for(i in parents){
    parname<-c(parname,paste0(i,"_",seq(1:ploidy)))
  }
  gen[,parname]<-NA

  #to obtain a vector with all homologue names.
  homol<-names(gen)
  homol<-homol[2:length(homol)]

  #for each marker, in each homologue a sampling procedure decides between A or B
  for(i in 1:length(markers[,1])){
    gen[i,-1]<-sample(c("A","B"),length(homol),replace=T)
  }

  #write it in a file using makeGen function.
  makeGen(output,gen)
}

#' Parental sampling
#' @description Generates a chromosome file based on sampling from a set of files containing parental
#' chromosomes (.dat files of PedigreeSim). To specify number of ancestrals, one can chose to
#' specify the number of ancestral groups, from how many different files should the parents be sampled.
#' By default, all files given will be used.
#'
#' @param files table files containing parental chromosomes, labelled as Parent_1,
#' Parent_2 ... Parent_Ploidy. Each file will constitute an "ancestral group"
#' @param ancgroups number of ancestral groups to use in the sampling. By default, all given ancestral groups
#' will be used.
#' @param parents number of parents to sample.
#' @param ploidy ploidy of the parents.
#'
#' @return
#' @export
#'
#' @examples
sampleChrom<-function(#Give a vector of file names, containing parental chromosomes. Returns a matrix
  #It generates only the most balanced situations (10 parents from 3 pops will give 3+3+4)
  files,
  ancgroups=length(files),
  parents=10,
  ploidy=4){

  repeat{
    chosen<-sample(files,ancgroups)
    if(length(unique(chosen))==ancgroups){break}
  }
  chosen<-sort(c(rep(chosen,parents%/%ancgroups),head(chosen,parents%%ancgroups)))

  #candidate parents of the chosen ancestral pools
  candparents<-sapply(unique(chosen),function(file){
    names<-as.matrix(data.table::fread(file,nrows=1,header=F))[-1]
    names<-unique(substr(names,1,nchar(names)-2))
  })

  #choose parents of each ancestral pool (without repeating)
  repeat{
    chosenparents<-apply(candparents[,chosen],2,sample,1)
    if(length(unique(chosenparents))==parents){break}
  }

  #Now we use the chosenparents vector to read the corresponding chromosome columns and paste them into one file
  result<-lapply(1:length(chosenparents),function(i){
    chroms<-paste0(chosenparents[i],"_",1:ploidy)
    data.table::fread(chosen[i],select=chroms,header=T)
  })
  result<-do.call(cbind,result)
  markers<-data.table::fread(chosen[1],select="marker",header=T)
  return(cbind(markers,result))
}


# Running the program ---------------

#' Run PedigreeSim from R studio
#' @description Run PedigreeSim from R studio by specifying a parameter file and the adress of
#' PedigreeSim.
#'
#' @param parfile name of the parameter file.
#' @param path name of the PedigreeSim program.
#'
#' @return
#' @export
#'
#' @examples
run.PedigreeSim <- function(
  parfile, #name of the parfile to be used
  path="PedigreeSIM/PedigreeSim.jar" #specify location of PedigreeSim file
) {

  #send the command to the cmd
  system2(command = "java",
          args = c("-jar",
                   path,
                   parfile))
}


# HWE populations ----------------------
hweGen <- function(
  map, #name of the map file to be used (with extension)
  parents, #number of parents
  name_parents = NULL,
  output=NULL, #name of the output
  ploidy=4,
  maf=NULL,
  object_map = F
){
  #read map file to obtain marker names.
  if(object_map){
    markers <- map
  }else{
    markers <- read.table(map,sep=" ",header = T)
  }


  #This distribution is made so that most
  #gene frequencies are between 0 and 0.2
  #to see underlying distribution do
  # quants <- seq(0,0.5,length.out = 100)
  # distro <- dbeta(quants,1.3,15)
  # plot(quants,distro/max(distro),type="l",
  #      main="MAF distribution (beta a=1.3, b=15)",
  #      xlab="MAF",ylab="frequency")
  if(is.null(maf)){
    maf =  rbeta(nrow(markers),1.3,15)
    maf[maf > 0.5] <- 0.5
  }else{
    if(!length(maf)==nrow(markers)){
      stop(length(maf)," marker frequencies provided for",
           nrow(markers)," markers in the map")
    }
  }

  gen <- t(sapply(maf,function(f){
    copies_A <- rbinom(parents*ploidy,size = 1,f) + 1
    return(c("A","B")[copies_A])
  }))

  #Create name vector for parents
  if(is.null(name_parents)){
    name_parents = paste0("P",0:(parents-1))
  }
  #Add suffixes for each chromosome
  name_parents <- sapply(name_parents, function(p) paste0(p,"_",1:ploidy))
  colnames(gen) <- name_parents

  #write it in a file using makeGen-style function.
  gen <- cbind(marker=as.character(markers[,1]),gen)
  write.table(gen,paste0(output,".gen"),quote=F,row.names = F)
}


hwe_plot <- function(ploidy, legend=T){

  #How many types of genotypes are there
  #i.e. for diploids 1,2,1 -> 1p^2 2pq 1q^2
  pascal <- choose(ploidy, 0:ploidy)
  plot(1,type="n",xlim=c(0,1),ylim=c(0,1),
       xlab= "MAF (freq of B)",ylab="Genotype frequency",
       main="HWE plot")
  cols <- viridis::viridis(ploidy+1)
  for(i in 0:ploidy){
    curve(pascal[i+1]*x^i*(1-x)^(ploidy-i),add=T,
          col=cols[i+1],lwd=2)
  }

  if(legend){
    leg <- sapply(0:ploidy,function(i){
      geno <- c("A","B")[c(rep(1,ploidy-i),rep(2,i))]
      paste(geno,collapse="")
    })
    legend("top",legend=leg,text.col=cols,bty="n")
  }

}

hwe_test <- function(
  ploidy = 4,
  mark_num = 1000,
  inds = 100
){

  #First we generate some allele frequencies,
  #in this case using a beta distribution
  maf =  rbeta(mark_num,1.3,15)
  maf[maf > 0.5] <- 0.5

  #We use the frequencies to obtain genotypes per chromosome
  gen <- t(sapply(maf,function(f){
    copies_A <- rbinom(inds*ploidy,size = 1,f) + 1
    return(c("A","B")[copies_A])
  }))

  #We group the chromosomes into single genotypes
  chrom_matrix <- sapply(1:inds,function(parent){
    chroms <- ((parent-1)*4+1):(parent*4)
    apply(gen[,chroms],1,function(alleles) paste(sort(alleles),collapse="") )
  })

  #We create genotype classes and check their frequency
  geno_classes <- sapply(0:ploidy,function(i){
    geno <- c("A","B")[c(rep(1,ploidy-i),rep(2,i))]
    paste(geno,collapse="")
  })
  freq_matrix <- apply(chrom_matrix,1,function(marker){
    sapply(geno_classes,function(g) sum(g==marker)/length(marker))
  })

  #We create a theoretical HWE plot
  hwe_plot(ploidy)
  #And we finish by plotting all the points
  cols <- viridis::viridis(ploidy+1)
  points(rep(maf,each=ploidy+1),as.vector(freq_matrix),
         col=rep(cols,mark_num),pch=19,cex=0.7)
}
