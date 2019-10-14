################################
### PedigreeSIM populations ####
################################
# With this file the ancestral populations and NAM populations can be simulated
source("mpQTL_fun.R")

#### Ancestral Simulations ####
#first, we define the adress that we want the Populations to be, in our case 10 pops
pops<-paste0("PedigreeSIM/Parents/Pop_",1:10)

founders<-10
popsize<-10
ngener<-49 #number of generations not including founders


#Make input files
for(name in pops){
  #Parameter files, indicating Potato chromosome and map files
  makePar(name,
          mapgen="PedigreeSIM/Potato", #mapfile to use
          chrom="PedigreeSIM/Potato") #chromosome file
  
  #Random pedigree files. Thanks to Rouland for his code in here
  makeRanPed(name,
             nfounder=founders, #number of founders
             popsize=popsize, #size of each generation
             ngener=ngener) #number of generations (not including founders)
  
  #Random genotypes for each of the founders
  ranGen(map="PedigreeSIM/Potato.map", #map to use as scheme
         parents=paste0("G0",1:10), #names of parents to simulate
         output=name, 
         ploidy=4)
}

#Run PedigreeSIM. PedSim messages will be printed in the R console
for(name in pops) run.PedigreeSim(paste0(name,".par"))

#### Parental Simulations ####
#First, we will reduce the genotypes files produced by Pedsim, as they contain ALL generations,
#and we only want the last one
parentgenotypes<-as.character(data.table::fread(paste0(pops[1],"_genotypes.dat"),nrows=1,header=F))

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

#Now we will sample the new files to make the Pedsim files of the NAM crosses
parentfiles<-paste0(pops,"_parents.dat") #names of parent genotype files
crosses<-100 #number of crosses of each NAM population
tablefile<-"PedigreeSIM/NAM_crosses/Table.txt"

#This couple lines create a table file that contains information about each cross
colnames<-c("cross_name","ancestors",paste0("parent_",0:9))
cat(c(colnames,"\n"),sep=" ",file = tablefile)

#In here we will make a set of folders, each of them with the cross files of 
#different NAM populations.
for(anc in 1:length(parentfiles)){
  #create a folder where to store the files
  folder<-paste0("PedigreeSIM/NAM_crosses/",anc,"_ancestral")
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
            mapgen="PedigreeSIM/Potato",
            chrom="PedigreeSIM/Potato",
            pedfile = file.path(folder,recordline[1]),
            founder = file.path(folder,recordline[1]))
  }
  
}

#Make a vector with the adresses of the crossfiles.par
crossfiles<-paste0("PedigreeSIM/NAM_crosses/",1:10,"_ancestral")
crossfiles<-sapply(1:length(crossfiles),function(file){
  paste0(crossfiles[file],"/cross",sprintf("%03d",0:(crosses-1)+crosses*(file-1)),".par")
  
})
crossfiles<-as.vector(crossfiles)

#run each file on PedigreeSim
for(file in crossfiles){
  run.PedigreeSim(file,
                  path="PedigreeSIM/PedigreeSim.jar")
}



