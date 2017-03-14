
######Co-occurrence pairwise correlations###### 

#Read in microbe data. Input datasets are datBacr3fotu3, datEukN5fotu3, datEukS4fotu3, datITS3fotu3. they were relativized prior to doubletons/singletons/summed<.2% removal. 
comm.dataEukS<-datEukS4fotu3
comm.dataEukN<-datEukN5fotu3
comm.dataBac<-datBacr3fotu3
comm.dataITS<-datITS3fotu3


#Read in plant data
plantcomp<-read.csv("/Users/farrer/Dropbox/Niwot Moving Uphill/Analysis/Niwot_MovingUpHill_comp2015.csv")
head(plantcomp)
names(plantcomp)[1]<-"Sample_name"

#Remove plant species only present in one or two plots; there are some plots that have plant data but not microbe data. 69 70 71 77 81 108 117 118 147 148 149 151. This is because when we started doing the surveys we were going to all plots for plants and only sample some for microbes, then we realized that that was insane!
dim(plantcomp)
plantcomp2<-plantcomp[,colSums(plantcomp>0)>2]
plantcomp2$LICHEN<-NULL
labelsPlant<-as.data.frame(cbind(otu=colnames(plantcomp2)[2:55],labels="Plant"))



#Merge things. all microbe datasets (not plants) should have the same 90 samples
#first merge comm.dataEuk with comm.data16S
#I need to remove all the description columns in one of the files, then merge
comm.dataEukSa<-cbind(Sample_name=comm.dataEukS$Sample_name,comm.dataEukS[,-c(1:31)])
comm.dataALL1<-merge(comm.dataBac,comm.dataEukSa,"Sample_name",sort=F,all.y=F,all.x=F)

comm.dataEukNa<-cbind(Sample_name=comm.dataEukN$Sample_name,comm.dataEukN[,-c(1:31)])
comm.dataALL2<-merge(comm.dataALL1,comm.dataEukNa,"Sample_name",sort=F,all.y=F,all.x=F)

comm.dataITSa<-cbind(Sample_name=comm.dataITS$Sample_name,comm.dataITS[,-c(1:31)])
comm.dataALL3<-merge(comm.dataALL2,comm.dataITSa,"Sample_name",sort=F,all.y=F,all.x=F)

comm.dataALL3$Sample_name


#then merge plants with microbes
comm.dataALL4<-merge(comm.dataALL3,plantcomp2,"Sample_name",sort=F,all.y=F)
comm.dataALL4$Sample_name
comm.dataALL<-comm.dataALL4[order(comm.dataALL4$Sample_name),]

#Make explanatory variable
comm.dataALL$lomehi<-factor(comm.dataALL$lomehi,levels=c("lo","me","hi"))
lomehiALL<-comm.dataALL$lomehi#30 lo, 29 me, 31 hi 
length(which(lomehiALL=="lo"))
length(which(lomehiALL=="me"))
length(which(lomehiALL=="hi"))

trts<-as.vector(levels(lomehiALL))

#see what plant diversity (richness) levels are in the lo me hi categories
comm.dataALL$Plant_Dens[which(comm.dataALL$lomehi=="lo")]# range 0-28
comm.dataALL$Plant_Dens[which(comm.dataALL$lomehi=="me")]# range 31-80
comm.dataALL$Plant_Dens[which(comm.dataALL$lomehi=="hi")]# range 81-739
comm.dataALL$Plant_Div[which(comm.dataALL$lomehi=="lo")]# range 0-8, mean 2.13
comm.dataALL$Plant_Div[which(comm.dataALL$lomehi=="me")]# range 3-14, mean 6.52
comm.dataALL$Plant_Div[which(comm.dataALL$lomehi=="hi")]# range 7-26, mean 14.90

cor.test(comm.dataALL$Plant_Dens,comm.dataALL$Plant_Div)


#Setup parallel backend to use 4 processors
cl<-makeCluster(4)
registerDoParallel(cl)
#start time
strt<-Sys.time()
results<-matrix(nrow=0,ncol=9)
options(warnings=-1)

for(a in 1:length(trts)){
  #pull the first element from the vector of treatments
  trt.temp<-trts[a]
  #subset the dataset for those treatments
  temp1<-subset(comm.dataALL, lomehi==trt.temp)####change to comm.data if not doing plants

  #to make it faster, I should take out zeros here. I am also taking out doubletons and singletons since i will never use them for their correlations and won't use them for correcting qvalue either
  temp<-cbind(temp1[,1:31],temp1[,((which(colSums(temp1[,32:dim(temp1)[2]]>0)>2))+31)])
  
  #in this case the community data started at column 32, so the loop for co-occurrence has to start at that point
  for(b in 32:(dim(temp)[2]-1)){
    
    results1<-foreach(c=(b+1):(dim(temp)[2]),.combine=rbind) %dopar% {
      species1.ab<-sum(temp[,b])
      species2.ab<-sum(temp[,c])
      species1.abfreq<-sum(temp[,b]>0)
      species2.abfreq<-sum(temp[,c]>0)

      if(species1.abfreq>0&species2.abfreq>0){#Ryan had this at 0
        test<-cor.test(temp[,b],temp[,c],method="spearman",na.action=na.rm)
        spearmanrho<-test$estimate
        spearmanp.value<-test$p.value
      }    
      if(species1.abfreq <=0 | species2.abfreq <= 0){
        spearmanrho<-0
        spearmanp.value<-1
      }	
      new.row<-c(trts[a],names(temp)[b],names(temp)[c],spearmanrho,spearmanp.value,species1.ab,species2.ab,species1.abfreq,species2.abfreq)
    }
    results<-rbind(results,results1)
  }
  print(a/length(trts))
}
head(results)
resultsold<-results
print(Sys.time()-strt) #Euk, bact, fungi, plants, took 1.3 days (31 hrs) with 4 cores
stopCluster(cl)

rownames(resultsold)<-1:dim(resultsold)[1]
results<-data.frame(resultsold[,1:3],(apply(resultsold[,4:dim(resultsold)[2]],2,as.numeric)))
names(results)<-c("trt","taxa1","taxa2","spearmanrho","spearmanp.value","ab1","ab2","ab1freq","ab2freq")
dim(results)

#from run with spearmanm, with a 3 sample frequency cutoff and removal of otus with less than .2% relative abundance
resultsS<-results



###### Edgelist creation #####

#I think the most correct way of doing this is: remove frequencies under whatever limit I want to set (10), do the qvalue correction on all spearman p values (positive and negative rhos), then subset and use only positive rho

edge_listsBEPa<-cbind(pairs=paste(resultsS$taxa1,resultsS$taxa2),resultsS) #this takes about 10 min
edge_listsBEPb<-subset(edge_listsBEPa,ab1freq>10&ab2freq>10) #this takes about 3 min
dim(edge_listsBEPb)
edge_listsBEPb$qval<-p.adjust(edge_listsBEPb$spearmanp.value,method="fdr")
edge_listsBEP<-subset(edge_listsBEPb,spearmanrho>0)
head(edge_listsBEP)
dim(edge_listsBEP)
min(edge_listsBEP$ab1freq)
dim(subset(edge_listsBEP,trt=="me"))

#how many plants and metazoa are in the dataset after different frequency cutoffs
listm<-c(as.character(edge_listsBEP$taxa1),as.character(edge_listsBEP$taxa2))
listm2<-unique(subset(listm,listm%in%labelsEukN4$otu)) #this doesn't list it twice if it is present in low and med
listp2<-unique(subset(listm,listm%in%labelsPlant2$otu)) #this doesn't list it twice if it is 
listm2
listp2

#with 9, you get more metazoa showing up in med (1) and low (1) but pattern is not nice
#with 10, you get fewer metazoa, but pattern is better
#with 11, the pattern is best but it is hard to explain an 11 cutoff

#see how many plant interactions are in the final network dataset
edge_listsBEPp<-edge_listsBEP[which(edge_listsBEP$taxa1%in%labelsEukN4$otu|edge_listsBEP$taxa2%in%labelsEukN4$otu),]
dim(edge_listsBEPm)
edge_listsBEPp<-edge_listsBEP[which(edge_listsBEP$taxa1%in%labelsPlant2$otu|edge_listsBEP$taxa2%in%labelsPlant2$otu),]
dim(edge_listsBEPp)
head(edge_listsBEPp)
temp<-subset(edge_listsBEPp,qval<.01&spearmanrho>.6&trt=="me")[,3:4]
temp
temp<-subset(edge_listsBEPp,qval<.01&spearmanrho>.6&trt=="hi")[,3:4]
temp
tax_table(datEukS4)[rownames(tax_table(datEukS4))=="denovo55021"]

#extract plant-microbe interactions, for Katie
edge_listsPlants1<-subset(edge_listsBEPa,ab1freq>5&ab2freq>5) 
ind<-labelsall[which(labelsall$labels=="Plant"),]
edge_listsPlants2<-edge_listsPlants1[which(edge_listsPlants1$taxa1%in%as.character(ind$otu)|edge_listsPlants1$taxa2%in%as.character(ind$otu)),]
dim(edge_listsPlants1)
dim(edge_listsPlants2)

edge_listsPlants2$qval<-p.adjust(edge_listsPlants2$spearmanp.value,method="fdr")
edge_listsPlants<-subset(edge_listsPlants2,spearmanrho>0)

