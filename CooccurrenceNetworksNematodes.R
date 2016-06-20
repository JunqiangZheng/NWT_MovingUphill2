#correlations including the 99% similarity nematode data (floating from 20g subsample and then sequenced)


#metazoa data: dat99Met2fotu2

comm.data99Met<-dat99Met2fotu2

#take out all sample data except Sample_name
comm.data99Meta<-cbind(Sample_name=comm.data99Met$Sample_name,comm.data99Met[,-c(1:31)])

#use comm.dataALL from the CooccurrenceNetworksEuksbac.R script and merge with nematode community data. There are 90 samples in common between the nematode and the bac/euk sequencing dataset
comm.dataALLm<-merge(comm.dataALL,comm.data99Meta,"Sample_name",sort=F,all.y=F,all.x=F)
comm.dataALLm$Sample_name
head(comm.dataALLm)[,1:34]

#same as before
trts


#only nematode-other interacitons, no nematode-nematode interactions
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
  temp1<-subset(comm.dataALLm, lomehi==trt.temp)####change to comm.data if not doing plants
  
  #to make it faster, I should take out zeros here. I am also taking out doubletons and singletons since i will never use them for their correlations and won't use them for correcting qvalue either
  temp2<-cbind(temp1[,1:31],temp1[,((which(colSums(temp1[,32:dim(temp1)[2]]>0)>2))+31)])
  ind<-which(substr(colnames(temp2),start=1,stop=2)=="md")
  tempmet<-temp2[,ind]
  temp<-temp2[,-ind]
  
  #take nematode data, for each nematode otu do correlations with every bacteria/plant otu. remember bact community data started at column 32, so the loop for co-occurrence has to start at that point
  for(b in 1:(dim(tempmet)[2])){
    results1<-foreach(c=32:(dim(temp)[2]),.combine=rbind) %dopar% {
      species1.ab<-sum(tempmet[,b])
      species2.ab<-sum(temp[,c])
      species1.abfreq<-sum(tempmet[,b]>0)
      species2.abfreq<-sum(temp[,c]>0)
      #I changed this so that it will calculate the correlation for all species pairs (except doubletons and singletons), and I can subset them later if I want to remove infrequent otus for qvalue calculation
      if(species1.abfreq>0&species2.abfreq>0){#Ryan had this at 0
        test<-cor.test(tempmet[,b],temp[,c],method="spearman",na.action=na.rm)
        spearmanrho<-test$estimate
        spearmanp.value<-test$p.value
      }    
      if(species1.abfreq <=0 | species2.abfreq <= 0){
        spearmanrho<-0
        spearmanp.value<-1
      }	
      new.row<-c(trts[a],names(tempmet)[b],names(temp)[c],spearmanrho,spearmanp.value,species1.ab,species2.ab,species1.abfreq,species2.abfreq)
      new.row		
    }
    results<-rbind(results,results1)
  }
  print(a/length(trts))
}
head(results)
resultsold<-results
print(Sys.time()-strt) #nematodes vs Euk, bact, plants, took 7.5min with 4 cores
stopCluster(cl)

rownames(resultsold)<-1:dim(resultsold)[1]
results<-data.frame(resultsold[,1:3],(apply(resultsold[,4:dim(resultsold)[2]],2,as.numeric)))
names(results)<-c("trt","taxa1","taxa2","spearmanrho","spearmanp.value","ab1","ab2","ab1freq","ab2freq")


head(results)

resultsM<-results

#####Nematode correlations######
edge_listsM<-cbind(pairs=paste(resultsM$taxa1,resultsM$taxa2),resultsM)
#edge_listsMb<-subset(edge_listsM,ab1freq>6&ab2freq>6)
#edge_listsMb$qval<-p.adjust(edge_listsMb$spearmanp.value,method="fdr")
#edge_listsMc<-subset(edge_listsMb,spearmanrho>0)
dim(edge_listsMc)

#subset metazoa 97% out from the euks list
edge_listsMf<-edge_listsM[-which(edge_listsM$taxa2%in%c(rownames(labelsEukN2))),]
head(edge_listsMf)
dim(edge_listsMf)
edge_listsMfb<-subset(edge_listsMf,ab1freq>5&ab2freq>5)
edge_listsMfb$qval<-p.adjust(edge_listsMfb$spearmanp.value,method="fdr")
edge_listsMfc<-subset(edge_listsMfb,spearmanrho>0)








#Including nematode-nematode interactions
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
  temp1<-subset(comm.dataALLm, lomehi==trt.temp)####change to comm.data if not doing plants
  
  #to make it faster, I should take out zeros here. I am also taking out doubletons and singletons since i will never use them for their correlations and won't use them for correcting qvalue either
  temp2<-cbind(temp1[,1:31],temp1[,((which(colSums(temp1[,32:dim(temp1)[2]]>0)>2))+31)])
  ind<-which(substr(colnames(temp2),start=1,stop=2)=="md")
  tempmet<-temp2[,ind]
  temp<-temp2[,-ind]
  
  resultsa<-matrix(nrow=0,ncol=9)
  resultsb<-matrix(nrow=0,ncol=9)
  
  #first correlate nematodes with nematodes
  for(b in 1:(dim(tempmet)[2]-1)){
    
    results1<-foreach(c=(b+1):(dim(tempmet)[2]),.combine=rbind) %dopar% {
      species1.ab<-sum(tempmet[,b])
      species2.ab<-sum(tempmet[,c])
      species1.abfreq<-sum(tempmet[,b]>0)
      species2.abfreq<-sum(tempmet[,c]>0)
      
      if(species1.abfreq>0&species2.abfreq>0){#Ryan had this at 0
        test<-cor.test(tempmet[,b],tempmet[,c],method="spearman",na.action=na.rm)
        spearmanrho<-test$estimate
        spearmanp.value<-test$p.value
      }    
      if(species1.abfreq <=0 | species2.abfreq <= 0){
        spearmanrho<-0
        spearmanp.value<-1
      }	
      new.row<-c(trts[a],names(tempmet)[b],names(tempmet)[c],spearmanrho,spearmanp.value,species1.ab,species2.ab,species1.abfreq,species2.abfreq)
    }
    resultsa<-rbind(resultsa,results1)
  }
  
  #take nematode data, for each nematode otu do correlations with every bacteria/plant otu. remember bact community data started at column 32, so the loop for co-occurrence has to start at that point
  for(b in 1:(dim(tempmet)[2])){
    results1<-foreach(c=32:(dim(temp)[2]),.combine=rbind) %dopar% {
      species1.ab<-sum(tempmet[,b])
      species2.ab<-sum(temp[,c])
      species1.abfreq<-sum(tempmet[,b]>0)
      species2.abfreq<-sum(temp[,c]>0)
      #I changed this so that it will calculate the correlation for all species pairs (except doubletons and singletons), and I can subset them later if I want to remove infrequent otus for qvalue calculation
      if(species1.abfreq>0&species2.abfreq>0){#Ryan had this at 0
        test<-cor.test(tempmet[,b],temp[,c],method="spearman",na.action=na.rm)
        spearmanrho<-test$estimate
        spearmanp.value<-test$p.value
      }    
      if(species1.abfreq <=0 | species2.abfreq <= 0){
        spearmanrho<-0
        spearmanp.value<-1
      }	
      new.row<-c(trts[a],names(tempmet)[b],names(temp)[c],spearmanrho,spearmanp.value,species1.ab,species2.ab,species1.abfreq,species2.abfreq)
      new.row		
    }
    resultsb<-rbind(resultsb,results1)
  }
  results<-rbind(results,resultsa,resultsb)
  print(a/length(trts))
}
head(results)
resultsold<-results
print(Sys.time()-strt) #nematodes vs Euk, bact, plants, took 7.7min with 4 cores
stopCluster(cl)

rownames(resultsold)<-1:dim(resultsold)[1]
results<-data.frame(resultsold[,1:3],(apply(resultsold[,4:dim(resultsold)[2]],2,as.numeric)))
names(results)<-c("trt","taxa1","taxa2","spearmanrho","spearmanp.value","ab1","ab2","ab1freq","ab2freq")


resultsM2<-results

##### Nematode correlations 2, including nematode-nematode interactions ######
edge_listsM2<-cbind(pairs=paste(resultsM2$taxa1,resultsM2$taxa2),resultsM2)

#take out metazoa 97% from the euks list
edge_listsM2a<-edge_listsM2[-which(edge_listsM2$taxa2%in%c(rownames(labelsEukN2))),]
head(edge_listsM2a)
dim(edge_listsM2a)
edge_listsM2b<-subset(edge_listsM2a,ab1freq>4&ab2freq>4)
edge_listsM2b$qval<-p.adjust(edge_listsM2b$spearmanp.value,method="fdr")
edge_listsM2c<-subset(edge_listsM2b,spearmanrho>0)

#take out nematodes 97% from the euks list and take out all metazoa except nematodes from the Met99 list. this edgelist is much smaller than above since we are taking taxa out of the 99% portion of the data (each of the 99% gets paired with ~4000 of the 97% taxa, so taking one 99% taxon out removes a lot of interactions)
ind<-labelsEukN2[which(labelsEukN2=="Nematoda")]
edge_listsM2d<-edge_listsM2[-which(edge_listsM2$taxa2%in%c(rownames(ind))),] #they are only in taxa 2

ind<-labels99Met2[which(labels99Met2!="Nematoda")]
edge_listsM2e<-edge_listsM2d[-which(edge_listsM2d$taxa2%in%c(rownames(ind))|edge_listsM2d$taxa1%in%c(rownames(ind))),] #they are in both taxa 2 and taxa 1

edge_listsM2f<-subset(edge_listsM2e,ab1freq>4&ab2freq>4)
edge_listsM2f$qval<-p.adjust(edge_listsM2f$spearmanp.value,method="fdr")
edge_listsM2g<-subset(edge_listsM2f,spearmanrho>0)

temp<-edge_listsM2g %>%
  filter(taxa1=="mdenovo226304"|taxa2=="mdenovo226304")
sort(temp$qval)    




