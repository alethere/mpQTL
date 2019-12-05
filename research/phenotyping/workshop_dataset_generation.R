#Dataset simulation for the workshop / mpQTL package
#' the purpose of this script is to generate a set of genotypic and phenotypic
#' data that can be used in the workshop of December 2019

#1 Create a reduced map, chromosomes 1 to 4 and unmapped markers ------------------------

map <- read.table("research/PedigreeSim/Potato.map",stringsAsFactors = F,header = T)
map <- map[map$chromosome < 5, ]
write.table(map,loc("workshop.map"),quote = F,row.names = F,col.names = T)

location <- "research/PedigreeSim/Workshop_data/"
loc <- function(text) paste0(location,text)

#2 Create a population that satisfies the segregation requirements of Yanlin ------------------
source("R/pedsim_fun.R")

makedir(location)
#1 First we create a pedigree for our dataset
founders <- paste0("founder_",1:20)

random_mixing <- function(founders,n_gen,pop_size,prefix=NULL,shift=NULL){
  #Function to create a dataframe of a "random mixing" population
  lineage <- list(data.frame(name = founders,
                             parent1 = NA,
                             parent2 = NA,
                             stringsAsFactors = F))
  ascii <- c(LETTERS,letters,0:9)

  get_code <- function(n) sapply(1:n,function(x)
    paste0(sample(ascii,8,replace=T),collapse = ""))
  get_parents <- function(n){
    sample(lineage[[i-1]][,1],n,replace=T)
  }

  for(i in 2:n_gen){
    parent1 <- get_parents(pop_size)
    parent2 <- get_parents(pop_size)

    self <- parent1 == parent2
    while(sum(self)>0){
      parent2[self] <- get_parents(sum(self))
      self <- parent1 == parent2
    }

    if(!is.null(shift)) gen <- i + shift
    else gen <- i

    lineage[[i]] <- data.frame(name = paste0(prefix,"G",gen,"_",get_code(pop_size)),
                            parent1, parent2,stringsAsFactors = F)
  }
  res <- do.call(rbind,lineage)
  return(res)
}
gen_size <- 50
gen_num <- 50; gen_A3 <- 30

#We create two parallel populations (ped_1 and ped_2)
#and a third population (ped_3) derived from generation 20 of ped_2
set.seed(7)
ped_1 <- random_mixing(founders,gen_num,gen_size,prefix = "A1")
ped_2 <- random_mixing(founders,gen_num,gen_size,prefix = "A2")

split_gen <- paste0("A2G",gen_num-gen_A3,"_")
ped_3 <- random_mixing(ped_2$name[grep(split_gen,ped_2$name)][1:10],
                       gen_A3,gen_size,prefix = "A3",shift = gen_num-gen_A3)

ped_2 <- ped_2[!grepl("founder",ped_2$name),]
ped_3 <- ped_3[!grepl("A2G[0-9]",ped_3[,"name"]),]


pedigree <- rbind(ped_1,ped_2,ped_3)
write.table(pedigree,loc("ancestral_workshop.ped"),
            col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#2 Then we create some gen files
hweGen(map,parents = 20,
       name_parents = founders,
       ploidy = 4,
       output = loc("founders_hwe"),object_map = T)
ranGen(map,parents = founders,
       output = loc("founders_random"),
       object_map = T,ploidy = 4)

source("R/viz_fun.R")
gen_hwe <- read.table(loc("founders_hwe.gen"),header = T)
gen_ran <- read.table(loc("founders_random.gen"),header = T)
source("R/mapQTL_fun.R")
colnames(gen)
K_ran <- calc.K(t(gen_ran[,-1]),haplotypes = T,ploidy = 4)
K_HwE <- calc.K(t(gen_hwe[,-1]),haplotypes = T,ploidy = 4)
pcoa.plot(K_ran)
pcoa.plot(K_HwE)


#Create parameter file
makePar(loc("workshop"),
        ploidy = 4,
        founder = loc("founders_hwe"),
        pedfile = loc("ancestral_workshop"),
        outname = loc("ancestral_workshop"))

#Create the chromosome file
lengths <- sapply(split(map$position,map$chromosome),max)
cens <- lengths/2
makeChrom(loc("workshop"),
          1:4,
          cM = lengths, centros = cens,
          pref = c(1,0.66,0.33,0),quad = 0)

run.PedigreeSim(loc("ancestral_workshop.par"),"research/PedigreeSim/PedigreeSim.jar")

#Generating the crosses
parents <-as.character(data.table::fread(loc("ancestral_workshop_genotypes.dat"),
                                                nrows=1,header=F))
parentcols <- grep(paste0("G",gen_num),parents)
genotypes <- data.table::fread(loc("ancestral_workshop_genotypes.dat"),
                               select=c(1,parentcols))

#Check of parent distances -----------------
K <- calc.K(t(genotypes[,-1]),haplotypes=T,ploidy = 4)
pcoa.plot(K,comp=c(1,2))
d <- as.dist(1-K)
attr(d,"Labels") <- substr(attr(d,"Labels"),1,2)
plot(ape::bionjs(d),type = "unrooted")


write.table(genotypes,file=loc("parents_workshop.gen"),quote = F,row.names = F,col.names = T)
found_alleles <- data.table::fread(loc("ancestral_workshop_founderalleles.dat"),
                                   select=c(1,parentcols))
write.table(found_alleles,file=loc("totallele_workshop.dat"),quote = F,row.names = F,col.names = T)
parents <- colnames(K)

parents[grepl("A1G50",parents)]
parents[grepl("A2G50",parents)]

#Parents per ancestral group
par_anc <- lapply(1:3,function(i) parents[grepl(paste0("A",i,"G50"),parents)])
set.seed(7)
sel_par <- sapply(par_anc,sample,3)

# 2.1 Cross design, structure 1 ---------------
#With these lines we create a mating scheme where each parent is used twice
#and they are mated the most connectedly possible.
index <- rep(1:length(sel_par),2)
crosses <- sapply(1:length(sel_par),function(i){
  res <- sel_par[index[i]]
  res2  <- sel_par[index[i+3]]
  return(c(res,res2))
})

n_per_cross <- sample(c(100,75,75,50,50,30,30,30,10),9)

ind <- as.vector(sel_par)
pedigree <- lapply(1:ncol(crosses),function(i){
  n <- n_per_cross[i]
  parent1 <- rep(crosses[1,i],n)
  parent2 <- rep(crosses[2,i],n)
  inds <- paste0("C",i,"_",1:n)
  return(data.frame(name = inds,parent1,parent2,stringsAsFactors = F))
})
pedigree <- do.call(rbind,pedigree)
pars <- data.frame(name = as.vector(sel_par),parent1 = NA,parent2 = NA,stringsAsFactors = F)
pedigree <- rbind(pars,pedigree)

#Now this pedigree table contains all crosses we're gonna make in our population
write.table(pedigree,loc("cross_workshop.ped"),quote = F,col.names = T,row.names = F)

index_sel_par <- colnames(genotypes) %in% as.vector(sapply(sel_par,paste0,"_",1:4))
index_sel_par[1] <- T
write.table(genotypes[,..index_sel_par],loc("parents_cross_workshop.gen"), quote=F,col.names = T,row.names = F)

makePar(loc("cross_workshop"),
        chrom=loc("workshop"),
        mapgen=loc("workshop"),
        founder = loc("parents_cross_workshop"),
        ploidy = 4)

run.PedigreeSim(loc("cross_workshop.par"),"research/PedigreeSim/PedigreeSim.jar")

# 2.2 Cross design structure 2 ----------------

makeNAMPed(as.vector(sel_par),offspring = 50,loc("NAM_workshop_2"))

makePar(loc("NAM_workshop_2"),
        chrom=loc("workshop"),
        mapgen=loc("workshop"),
        pedfile = loc("NAM_workshop_2"),
        founder = loc("parents_cross_workshop"),
        ploidy = 4)

run.PedigreeSim(loc("NAM_workshop_2.par"),"research/PedigreeSim/PedigreeSim.jar")

# 2.3 Cross design structure 3 ---------------------

makeChrom(loc("Non_pref_workshop"),
          1:4,
          cM = lengths, centros = cens,
          pref = rep(0,4),quad = 0)

makeNAMPed(as.vector(sel_par),offspring = 50, loc("NAM_workshop_nopref"))

write.table(genotypes[,..index_sel_par],
            loc("parents_cross_workshop_nopref.gen"), quote=F,col.names = T,row.names = F)


makePar(loc("NAM_workshop_nopref"),
        chrom=loc("Non_pref_workshop"),
        mapgen=loc("workshop"),
        pedfile = loc("NAM_workshop_nopref"),
        founder = loc("parents_cross_workshop_nopref"),
        ploidy = 4)

run.PedigreeSim(loc("NAM_workshop_nopref.par"),"research/PedigreeSim/PedigreeSim.jar")



# 3 Cross analysis ---------------
geno <- read.table("research/PedigreeSim/Workshop_data/cross_workshop_alleledose.dat",header = T)
K <- calc.K(t(geno[,-1]))
pop <- substr(colnames(geno)[-1],1,2)
pcoa.plot(K,col=pop,h=c(0,360))

id <- sapply(unique(pop),function(p) p == pop)
pch <- id%*%matrix(1:ncol(id),ncol=1)

cols <- select.col(ncol(id),coltype = "sequential")[pch]
pc <- prcomp(K)
scatterplot3d::scatterplot3d(x = pc$rotation[,1], y = pc$rotation[,2], z = pc$rotation[,3],
                             xlab = var.pc[1],ylab=var.pc[2],
                             zlab=var.pc[3],pch = pch,
                             color = cols,type = "p")

plotly::plot_ly(x = pc$rotation[,1], y = pc$rotation[,2], z = pc$rotation[,3],
                name = pop,
                type = "scatter3d",
                mode = "marker")


# Phenotyping ---------------
source("R/pheno_fun.R")
found <- link_NAM(loc("cross_workshop_founderalleles.dat"),
                        totallele = loc("totallele_workshop.dat"),
                        parents = 9,
                        ploidy = 4)

allele_count <- apply(found,1,function(g) length(unique(g)) )
hist(allele_count,20) #We have a high number of alleles

set.seed(7)
QTLpos <- sort(sample(1:nrow(map),2))
map[QTLpos,]

#Take out some markers to put them in "unmapped"
set.seed(13)
unmapped <- sample((1:nrow(map))[-QTLpos],150)
chr_0 <- data.frame(marker = map$marker[unmapped],
                    chromosome = 0,
                    position = 1:length(unmapped))

new_map <- rbind(map[-unmapped,],chr_0)
new_geno <- rbind(geno[-unmapped,],geno[unmapped,])
new_found <- rbind(found[-unmapped,],found[unmapped,])

#' Creating phenotypes is complicated, because assigning gene effects
#' is quite arbitrary and can take many forms. If one looks at the prevalence of
#' different alleles accross the randomly chosen parents (see below, dosage.X),
#' you can see there are some alleles unique to ancestral group, and some
#' alleles shared between them. In order to assign some effects, we can take the
#' evolutionary perspective: there were some original wt alleles that conferred
#' an additive effect X.
#' In our simulation, A1 could represent a population where
#' the trait has been selected against, and thus the common alleles have lower effects.
#' In A2, there is a wide range of ecotypes, so the alleles have a wide range of effects.
#' In A3, a subsample of ecotypes from A2 was taken, and created a new lineage where
#' most alleles produce a high phenotype.
#' Thus, alleles common between ancestral groups must be some form of wt,
#' while exclusive alleles can follow the group behaviour.
#'

View(dosage.X(found[QTLpos[1],],haplotype = T,ploidy=4))
dosage.X(found[QTLpos[2],],haplotype = T,ploidy=4)

anc_alleles_1 <- lapply(1:3,function(i){
  g <- paste0("A",i)
  a <- grepl(g,colnames(found)[1:36])
  sort(unique(found[QTLpos[1],1:36][a]))
})
names(anc_alleles_1) <- paste0("A",1:3)

anc_alleles_2 <- lapply(1:3,function(i){
  g <- paste0("A",i)
  a <- grepl(g,colnames(found)[1:36])
  sort(unique(found[QTLpos[2],1:36][a]))
})
names(anc_alleles_2) <- paste0("A",1:3)

anc_alleles <- c(anc_alleles_1, anc_alleles_2)

anc_mean <- c(40,60,30,50,20,20,50)
names(anc_mean) <- c("1","2","3","1;2","2;3","1;3","1;2;3")
anc_sd <- c(20,5,30,2,10,2,2)/2
names(anc_sd) <- c("1","2","3","1;2","2;3","1;3","1;2;3")

effects <- lapply(1:2, function(i){
  Qgenes <- found[QTLpos[i],]
  anc_alleles <- lapply(1:3,function(i){
    g <- paste0("A",i)
    a <- grepl(g,colnames(found))
    sort(unique(Qgenes[a]))
  })
  names(anc_alleles) <- paste0("A",1:3)
  anc_matrix <- sapply(anc_alleles,function(a) unique(Qgenes) %in% a)
  grouping <- apply(anc_matrix,1,function(a) paste(which(a),collapse=";") )

  unlist(effect_gen(Qgenes,anc_alleles, seed = i+-2091,
                    anc_sd = anc_sd[sort(unique(grouping))],
                    anc_mean = anc_mean[sort(unique(grouping))],
                    method = "anc_average"
                    ))
})

herits <- c(0.2,0.4,0.6,0.8)
phe <- sapply(herits,function(h){
  pheno(found[QTLpos,],effects,ploidy = 4,herit = h,mu = 50,seed = 5)
})
colnames(phe) <- herits

cof <- rep(1:3,length.out = nrow(phe))



pop_greenhouse <- split(unique(pop),rep(1:3,length.out=length(unique(pop))))
cof_mat <- sapply(pop_greenhouse,function(p) pop %in% p)
means <- c(5,-10,15)
greenhouse_effect <- apply(cof_mat,1,function(x) rnorm(1,mean=means[which(x)]))

greenhouse_effect <- sapply(cof,function(x) rnorm(1,mean=means[x]))


boxplot(phe[,3] + greenhouse_effect~substr(names(phe[,2]),1,2))
phenotypes <- cbind(phe,phe + greenhouse_effect)
colnames(phenotypes) <- c(0.2,0.4,0.6,0.8,paste0("cof_",c(0.2,0.4,0.6,0.8)))
rownames(phenotypes) <- gsub("_1$","",rownames(phenotypes))
rownames(phenotypes) <- gsub("_$","",rownames(phenotypes))

cofa <- apply(cof_mat,1,which)

phe <- phe[,c(3,7)]
colnames(phe) <- paste0("phenotype",1:2)

#4 Data saving -------------------
res <- map.QTL(phe,genotypes =  geno[,-1],map = map,ploidy = 4,K = T)
data <- list(pheno = phe,
             map = map,
             dosage = geno[,-1],
             founder = found,
             cofactor = cofa,
             result = res)
saveRDS(data,"vignette/workshop_data.RDS")

new_res <- map.QTL(phe,genotypes =  new_geno[,-1],map = new_map,ploidy = 4,K = T)
new_data <- list(pheno = phe,
             map = new_map,
             snp = new_geno[,-1],
             founder = new_found,
             cofactor = cofa,
             result = new_res)
saveRDS(new_data,"vignette/new_workshop_data.RDS")

effect_gen <- function(
  gen, #genotypes (per chromosome of individual)
  anc_alleles, #matrix (per column)/list defining which alleles belong to which ancestral group
  anc_sd=1,
  seed=7,
  method = "anc_average", #or "fix_eff"
  size = 3, #size of the effect per ancestral group
  n = 1, #number of alleles with an effect
  anc_mean = NULL
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
    grouping <- apply(anc_matrix,1,function(b) paste(which(b),collapse=";"))
    grouped_a <- split(a,grouping)
    anc_n <- sapply(grouped_a,length)

      # Average per anc group routine -----
      #Define the ancestral averages
      if(method == "anc_average"){
        set.seed(seed*k)

        if(is.null(anc_mean))  anc_mean <- runif(length(anc_n),0,1)
        else if(length(anc_mean) != length(anc_n)) stop("Number of anc_means provided is not equal to number of anc_groups")

        if(length(anc_sd) < length(anc_n)) anc_sd <- rep(anc_sd,length(anc_n))

        eff <- lapply(1:length(anc_n),function(i){
          set.seed(seed*k*i)
          #randomly generate effects with a specific seed
          e <- rnorm(anc_n[i],mean = anc_mean[i], sd=anc_sd[i])
          #which group are we taking (to put the right allele names)
          group <- names(grouped_a)[i]
          names(e) <- names(grouping)[grouping == group]
          return(e)
        })
        names(eff) <- names(anc_n)

      }else if(method == "fix_eff"){
        #Specific size per allele ------
        if(length(anc_n) > length(size)) size <- rep(size,length(anc_n))
        eff <- lapply(1:length(anc_n), function(i){
          e <- c(rep(size[i],n),rep(0,anc_n[i]-n))

          #which group are we taking (to put the right allele names)
          group <- names(grouped_a)[i]
          names(e) <- names(grouping)[grouping == group]
          return(e)
        })
        names(eff) <- names(anc_n)
      }

      res <- do.call(c,eff)
      names(res) <- unlist(grouped_a)
      return(res)
    })

  return(effects)
}


# Pedigree Plotting ----------------

ped2net <- function(pedigree){
  parents <- unique(unlist(pedigree[,c(2,3)]))
  parents <- na.omit(parents)

  net <- lapply(parents,function(p){
    pedigree$name[(pedigree$parent1 %in% p| pedigree$parent2 %in% p)]
  })
  names(net) <- parents
}
