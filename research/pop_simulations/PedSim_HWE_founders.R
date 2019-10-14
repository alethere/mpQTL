

# PedigreeSim with realistic gene frequencies ---------------
# Alejandro Thérèse Navarro
# July 2019

#' In this file we try to recreate Ancestral Populations and NAM populations,
#' but starting the simulations with genotypes defined by HW equilibrium

#1 Functions -----------------
source("pedsim_fun.R")


#test HWE -------------
png("Plots/MAF distribution for HWE simulations.png",
    res=300, width = 2000, height = 2000)
quants <- seq(0,0.5,length.out = 100)
distro <- dbeta(quants,1.3,15)
plot(quants,distro/max(distro),type="l",
     main="MAF distribution (beta a=1.3, b=15)",
     xlab="MAF",ylab="frequency")
dev.off()
png("Plots/HWE_test_beta_distribution.png",
    res=250,width = 2000,height = 1400)
hwe_test(4,1000,300)
dev.off()


#2 Main -------------

### Ancestral groups #### 

# 2.1 Parameter definition -----------
founders<-10
popsize<-100
ngener<-49 #number of generations not including founders
ploidy <- 4
ancestral_pops <- 3 #number of ancestral pops
folder_name <- "PedigreeSim/HWE_common" #name of the folder where results will be stored
gen_base <- "PedigreeSim/Potato"
map <- paste0(gen_base,".map")

# 2.2 Folder definition --------------
parent_folder <- paste0(folder_name,"/Parents/")
NAM_folder <- paste0(folder_name,"/NAM/")
for(f in c(folder_name,parent_folder,NAM_folder)) makedir(f)
pops <- paste0(parent_folder,"Pop_",1:ancestral_pops)

marknum <- nrow(data.table::fread(map))
maf <- rbeta(marknum,1.3,15)
maf[maf > 0.5] <- 0.5

#Make input files
for(name in pops){
  #Parameter files, indicating Potato chromosome and map files
  makePar(name,
          mapgen = gen_base, #mapfile to use
          chrom = gen_base) #chromosome file
  
  #Random pedigree files. Thanks to Roland for his code in here
  makeRanPed(name,
             nfounder=founders, #number of founders
             popsize=popsize, #size of each generation
             ngener=ngener) #number of generations (not including founders)
  
  #Random genotypes for each of the founders
  hweGen(map=map, #map to use as scheme
         parents=founders, #names of parents to simulate
         name_parents = paste0("G0_",sprintf("%02d",1:founders)),
         output=name, 
         ploidy=4,
         maf = maf)
}
#Run PedigreeSIM. PedSim messages will be printed in the R console
for(name in pops) run.PedigreeSim(paste0(name,".par"))

#### Parental Simulations ####
#First, we will reduce the genotypes files produced by Pedsim, as they contain ALL generations,
#and we only want the last one
parentgenotypes<-as.character(data.table::fread(paste0(pops[1],"_genotypes.dat"),
                                                nrows=1,header=F))

#we will create a file with only the genotypes of the parents, and change the parents names
parentcols<-parentgenotypes[(2+founders*ploidy+(ngener-1)*popsize*ploidy): #start of parents
                              (1+founders*ploidy+ngener*popsize*ploidy)] #last column

#We now create the new file with new names for the parents
#We will call them A0P00 for the first parent of the first ancestral group.
for(name in pops){
  genotypes<-data.table::fread(paste0(name,"_genotypes.dat"),select=parentcols)
  markers<-data.table::fread(paste0(name,"_genotypes.dat"),select="marker")
  newnames<-paste0("A",(which(name==pops)-1),"P",sprintf("%02.f",0:(popsize-1)))
  newnames<-sort(as.vector(outer(newnames,paste0("_",1:4),FUN="paste0")))
  colnames(genotypes)<-c(newnames)
  
  genotypes<-cbind(markers,genotypes)
  write.table(genotypes,file=paste0(name,"_parents.dat"),quote = F,row.names = F)
}

#Make the Totallele file for the linkNAM function
ancestors <- list()
for(i in 1:length(pops)){
  name <- pops[i]
  genotypes<-data.table::fread(paste0(name,"_founderalleles.dat"),select=parentcols)
  newnames<-paste0("A",i-1,"P",sprintf("%02.f",0:(popsize-1)))
  newnames<-sort(as.vector(outer(newnames,paste0("_",1:4),FUN="paste0")))
  colnames(genotypes)<-c(newnames)
  
  ancestors[[name]] <- genotypes+ploidy*founders*(i-1)
}
markers<-data.table::fread(paste0(name,"_genotypes.dat"),select="marker")
names(ancestors) <- NULL
ancestors <- do.call(cbind,ancestors)
write.table(cbind(markers,ancestors),
            file = paste0(parent_folder,"Totallele.txt"),
            quote = F,col.names = T,row.names = F)

#Now we will sample the new files to make the Pedsim files of the NAM crosses
parentfiles<-paste0(pops,"_parents.dat") #names of parent genotype files
crosses <- 10 #number of crosses of each NAM population
tablefile <- paste0(NAM_folder,"Table.txt")

#This couple lines create a table file that contains information about each cross
colnames <-c("cross_name","ancestors",paste0("parent_",0:9))
cat(c(colnames,"\n"),sep=" ",file=tablefile)

#In here we will make a set of folders, each of them with the cross files of 
#different NAM populations.
for(anc in 1:length(parentfiles)){
  #create a folder where to store the files
  folder<-paste0(NAM_folder,anc,"_ancestral")
  makedir(folder)
  
  for(nam in 1:crosses){#for each cross
    recordline<-c()#the record line to be printed in the table
    count<-((anc-1)*crosses+nam)#the total number of iterations that we are at
    recordline[1]<-paste0("cross",sprintf("%03.f",count-1))#name of the cross
    recordline[2]<-anc #ancestral number
    
    #Select parents, always choosing the most balanced combination (i.e 3 ancestrals: 3,3,4)
    parents<-sampleChrom(parentfiles,anc,ploidy=4)#select parents
    
    #Read the names of the selected 9 parents
    recordline<-c(recordline,
                  unique(substr(names(parents)[-1],0,nchar(names(parents)[-1])-2)))
    
    cat(recordline,"\n",file=tablefile,append=T,sep=" ")#write into file
    
    #make pedigree file
    makeNAMPed(parents=recordline[-1:-2],
               offspring=50,
               output=file.path(folder,recordline[1]),
               prefix=paste0(sprintf("%03.f",count-1),"_")#prefix to add to the offspring
    )
    
    #Write the genetic file
    makeGen(data=parents,
            output=file.path(folder,recordline[1]))
    
    #Write the parameter file
    makePar(output=file.path(folder,recordline[1]),
            mapgen=gen_base,
            chrom=gen_base,
            pedfile = file.path(folder,recordline[1]),
            founder = file.path(folder,recordline[1]))
  }
  
}


#Make a vector with the adresses of the crossfiles.par
crossfiles<-paste0(NAM_folder,1:3,"_ancestral")
crossfiles<-sapply(1:length(crossfiles),function(file){
  paste0(crossfiles[file],"/cross",sprintf("%03d",0:(crosses-1)+crosses*(file-1)),".par")
  
})
crossfiles<-as.vector(crossfiles)

#run each file on PedigreeSim
for(file in crossfiles){
  run.PedigreeSim(file,
                  path="PedigreeSim/PedigreeSim.jar")
}