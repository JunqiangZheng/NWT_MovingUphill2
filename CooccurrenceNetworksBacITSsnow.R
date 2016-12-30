
### Cooccurrence Networks across a snowdepth gradient for NSF preproposal ###

save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill2_Workspace_Analysis_Snow.Rdata")  #

load("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill2_Workspace_Analysis_Snow.Rdata")  #


#I had a problem with Git, I did a commit, but then it wouldn't sync because a .RDataTmp file was over 100GB in size. I did the following in terminal to uncommit. In the future I need to not stage/commit that .RDataTmp file
git commit -m “trying to undo my last commit when has a huge file in it“
git reset HEAD~
  

#Read in microbe data. Input datasets are datBacr3fotu3, datEukN5fotu3, datEukS4fotu3, datITS3fotu3. they were relativized prior to doubletons/singletons/summed<.2% removal. 
#comm.dataEukS<-datEukS4fotu3
#comm.dataEukN<-datEukN5fotu3
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
comm.dataITSa<-cbind(Sample_name=comm.dataITS$Sample_name,comm.dataITS[,-c(1:31)])
comm.dataALL1<-merge(comm.dataBac,comm.dataITSa,"Sample_name",sort=F,all.y=F,all.x=F)

comm.dataALL1$Sample_name

#then merge plants with microbes
comm.dataALL4<-merge(comm.dataALL1,plantcomp2,"Sample_name",sort=F,all.y=F)
comm.dataALL4$Sample_name
comm.dataALL5<-comm.dataALL4[order(comm.dataALL4$Sample_name),]

#Make explanatory variable, bring in snowdepth data
mapsnow<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/ITS_Niwot_20072015_All_MapFilenewlomehisnow.txt")
mapsnow2<-data.frame(mapsnow[,c(1,31)])

comm.dataALL6<-base::merge(comm.dataALL5,mapsnow2,"X.SampleID")
comm.dataALL7<-data.frame(comm.dataALL6[,c(1:31)],snowdepth=comm.dataALL6[,dim(comm.dataALL6)[2]],comm.dataALL6[,32:(dim(comm.dataALL6)[2]-1)])
comm.dataALL7[1:10,1:33]

#take out plots where there are 0 plants
comm.dataALL8<-subset(comm.dataALL7,Plant_Dens>0)

#this classification is very different from the one in the 2007/2015 analysis becuase here there are many more plots with deep snow, so it got reshuffled differently
hist(comm.dataALL8$snowdepth)
snow<-ifelse(comm.dataALL8$snowdepth<145,"lo",NA) #26 lo
length(which(snow=="lo"))
ind<-which(comm.dataALL8$snowdepth>=145&comm.dataALL8$snowdepth<250) #26 me
length(ind)
snow[ind]<-"me"
ind<-which(comm.dataALL8$snowdepth>=250&comm.dataALL8$snowdepth<500) #29 hi
length(ind)
snow[ind]<-"hi"

comm.dataALL8$lomehi<-factor(snow,levels=c("lo","me","hi"))
comm.dataALL<-comm.dataALL8
comm.dataALL$Description<-NULL
lomehi<-comm.dataALL8$lomehi
lomehiALL<-comm.dataALL$lomehi
length(which(lomehiALL=="lo"))
length(which(lomehiALL=="me"))
length(which(lomehiALL=="hi"))

trts<-as.vector(levels(lomehiALL))

comm.dataALLplant<-cbind(lomehi=lomehiALL,comm.dataALL[,4929:4982])
comm.dataALLmicrobe<-cbind(lomehi=lomehiALL,comm.dataALL[,32:4928])




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
  tempp1<-subset(comm.dataALLplant, lomehi==trt.temp)
  tempm1<-subset(comm.dataALLmicrobe, lomehi==trt.temp)
  
  #to make it faster, I should take out zeros here. I am also taking out doubletons and singletons since i will never use them for their correlations and won't use them for correcting qvalue either
  tempp<-cbind(lomehi=tempp1[,1],tempp1[,((which(colSums(tempp1[,2:dim(tempp1)[2]]>0)>2))+1)])
  tempm<-cbind(lomehi=tempm1[,1],tempm1[,((which(colSums(tempm1[,2:dim(tempm1)[2]]>0)>2))+1)])
  
  #in this case the community data started at column 2, so the loop for co-occurrence has to start at that point
  for(b in 2:(dim(tempp)[2])){
    
    results1<-foreach(c=(2):(dim(tempm)[2]),.combine=rbind) %dopar% {
      species1.ab<-sum(tempp[,b])
      species2.ab<-sum(tempm[,c])
      species1.abfreq<-sum(tempp[,b]>0)
      species2.abfreq<-sum(tempm[,c]>0)
      
      if(species1.abfreq>0&species2.abfreq>0){#Ryan had this at 0
        test<-cor.test(tempp[,b],tempm[,c],method="spearman",na.action=na.rm)
        spearmanrho<-test$estimate
        spearmanp.value<-test$p.value
      }    
      if(species1.abfreq <=0 | species2.abfreq <= 0){
        spearmanrho<-0
        spearmanp.value<-1
      }	
      new.row<-c(trts[a],names(tempp)[b],names(tempm)[c],spearmanrho,spearmanp.value,species1.ab,species2.ab,species1.abfreq,species2.abfreq)
    }
    results<-rbind(results,results1)
  }
  print(a/length(trts))
}
head(results)
resultsold<-results
print(Sys.time()-strt) #took 2 minutes
stopCluster(cl)

rownames(resultsold)<-1:dim(resultsold)[1]
results<-data.frame(resultsold[,1:3],(apply(resultsold[,4:dim(resultsold)[2]],2,as.numeric)))
names(results)<-c("trt","taxa1","taxa2","spearmanrho","spearmanp.value","ab1","ab2","ab1freq","ab2freq")
dim(results)



#from run with spearmanm, with a 3 sample frequency cutoff and removal of otus with less than .2% relative abundance
resultsS<-results



###### Edgelist creation #####

#I think the most correct way of doing this is: remove frequencies under whatever limit I want to set (5), do the qvalue correction on all spearman p values (positive and negative rhos), then subset and use only positive rho

edge_listsBEPa<-cbind(pairs=paste(resultsS$taxa1,resultsS$taxa2),resultsS) #this takes about 10 min
edge_listsBEPb<-subset(edge_listsBEPa,ab1freq>5&ab2freq>5) #this takes about 3 min
dim(edge_listsBEPb)
edge_listsBEPb$qval<-p.adjust(edge_listsBEPb$spearmanp.value,method="fdr")
#edge_listsBEP<-subset(edge_listsBEPb,spearmanrho>0)

#take moss out
edge_listsBEPc<-subset(edge_listsBEPb,taxa1!="MOSS")

edge_listsBEP<-edge_listsBEPc
head(edge_listsBEP)
dim(edge_listsBEP)
min(edge_listsBEP$ab1freq)
dim(subset(edge_listsBEP,trt=="me"))








#####Plant-microbe interactions#####
#use edge_listsBEP, cutoff frequency is 5

#Low snow
inputlo<-subset(edge_listsBEP,qval<.05&abs(spearmanrho)>.5&trt=="lo")[,3:4]
dim(inputlo)
inputlo2<-subset(edge_listsBEP,qval<.05&abs(spearmanrho)>.5&trt=="lo")[,c(3:5)]
graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otu"
colorgraph3<-merge(verticesgraph3,labelsall,"otu",all.y=F,all.x=F,sort=F)
#change all bacteria colors to purple, and plants to green, fungi to pink?
colorgraph3$color[which(colorgraph3$group=="PhotosyntheticBacteria")]<-"#7879BC"
colorgraph3$color[which(colorgraph3$group=="Plant")]<-"#94BA3C"
sizegraph3<-ifelse(colorgraph3$labels=="Plant",8,6)#was 6,4
shapegraph3<-ifelse(colorgraph3$labels=="Plant","square","circle")
linetype3<-ifelse(inputlo2$spearmanrho>0,1,2)
linecolor3<-ifelse(inputlo2$spearmanrho>0,"black","#E95275")
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/networklowsnow.pdf")
plot(graph3,vertex.size=sizegraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=F,vertex.label=NA,vertex.shape=shapegraph3,edge.width=2,edge.color=linecolor3)#,edge.lty=linetype3,,edge.color="gray40"
dev.off()
sum(linetype3==1) #positive interactions
sum(linetype3==2) #negative interactions
sum(linetype3==2)/(sum(linetype3==1)+sum(linetype3==2))

#Medium snow
inputme<-subset(edge_listsBEP,qval<.05&abs(spearmanrho)>.5&trt=="me")[,3:4]
dim(inputme)
inputme2<-subset(edge_listsBEP,qval<.05&abs(spearmanrho)>.5&trt=="me")[,3:5]
graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otu"
colorgraph2<-merge(verticesgraph2,labelsall,"otu",all.y=F,all.x=F,sort=F)
#change all bacteria colors to purple, and plants to green, fungi to pink?
colorgraph2$color[which(colorgraph2$group=="PhotosyntheticBacteria")]<-"#7879BC"
colorgraph2$color[which(colorgraph2$group=="Plant")]<-"#94BA3C"
#colorgraph2$color[which(colorgraph2$group=="Fungi")]<-"#E95275"
sizegraph2<-ifelse(colorgraph2$labels=="Plant",8,6)#was 6,4
shapegraph2<-ifelse(colorgraph2$labels=="Plant","square","circle")
linetype2<-ifelse(inputme2$spearmanrho>0,1,2)
linecolor2<-ifelse(inputme2$spearmanrho>0,"black","#E95275")
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/networkmedsnow.pdf")
plot(graph2,vertex.size=sizegraph2,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=F,vertex.label=NA,vertex.shape=shapegraph2,edge.width=2,edge.color=linecolor2)#,edge.color="gray40",edge.lty=linetype2
dev.off()
sum(linetype2==1) #positive interactions
sum(linetype2==2) #negative interactions
sum(linetype2==2)/(sum(linetype2==1)+sum(linetype2==2))

#High density
inputhi<-subset(edge_listsBEP,qval<.05&abs(spearmanrho)>.5&trt=="hi")[,3:4]
dim(inputhi)
inputhi2<-subset(edge_listsBEP,qval<.05&abs(spearmanrho)>.5&trt=="hi")[,3:5]
graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otu"
colorgraph1<-merge(verticesgraph1,labelsall,"otu",all.y=F,all.x=F,sort=F)
#change all bacteria colors to purple, and plants to green, fungi to pink?
colorgraph1$color[which(colorgraph1$group=="PhotosyntheticBacteria")]<-"#7879BC"
colorgraph1$color[which(colorgraph1$group=="Plant")]<-"#94BA3C"
#colorgraph1$color[which(colorgraph1$group=="Fungi")]<-"#E95275"
sizegraph1<-ifelse(colorgraph1$labels=="Plant",8,6)#was 6,4
shapegraph1<-ifelse(colorgraph1$labels=="Plant","square","circle")
linetype1<-ifelse(inputhi2$spearmanrho>0,1,2)
linecolor1<-ifelse(inputhi2$spearmanrho>0,"black","#E95275")
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/networkhisnow.pdf")
plot(graph1,vertex.size=sizegraph1,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=F,vertex.shape=shapegraph1,edge.width=2,vertex.label=NA,edge.color=linecolor1)#,edge.lty=linetype1,edge.color="gray40"
dev.off()
sum(linetype1==1) #positive interactions
sum(linetype1==2) #negative interactions
sum(linetype1==2)/(sum(linetype1==1)+sum(linetype1==2))

#legend
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/networklegend.pdf")
plot(c(1,1),c(1,1))
legend("topleft",c("Plant","Bacteria","Fungi"),pt.bg=c("#94BA3C","#7879BC","#F6EC32"),pch=c(22,21,21),bty="n",cex=1.4)
legend("bottomleft",c("Positive interaction","Negative interaction"),lty=1,col=c(1,"#E95275"),lwd=2,bty="n",cex=1.4)
dev.off()


#potential interactions
testlo<-subset(edge_listsBEP,trt=="lo")
unique(testlo$taxa1)
dim(testlo)
testme<-subset(edge_listsBEP,trt=="me")
unique(testme$taxa1)
dim(testme)
testhi<-subset(edge_listsBEP,trt=="hi")
unique(testhi$taxa1)
dim(testhi)

#strenght of interaction
mean(abs(inputlo2$spearmanrho))
mean(abs(inputme2$spearmanrho))
mean(abs(inputhi2$spearmanrho))

mean(inputlo2$spearmanrho)
mean(inputme2$spearmanrho)
mean(inputhi2$spearmanrho)

mean(abs(testlo$spearmanrho))
mean(abs(testme$spearmanrho))
mean(abs(testhi$spearmanrho))
std.error(abs(testlo$spearmanrho))
std.error(abs(testme$spearmanrho))
std.error(abs(testhi$spearmanrho))






#### Pulling out DM and SB for networks ####

#Clustering analysis on the plant community data
#ward.D with k=8
#ward.D2
plantdist<-dist(comm.dataALLplant[,-1],method="manhattan")
plantclust<-hclust(plantdist,method="ward.D")#ward.D2
plot(plantclust,hang=-1)
groups <- cutree(plantclust, k=8)
groups
cbind(groups,comm.dataALLplant[,c("lomehi","GEUROS","SILACA","CARPYR","JUNDRU","TRIDAS","MINOBT")])
m1<-metaMDS(comm.dataALLplant[,-1], distance="bray", k=2, autotransform=F,trymax=200)
spps  <- scores(m1, display = "species")
plot(scores(m1),col=groups)
text(spps, labels = rownames(spps), col = "black", cex = 0.6)#ifelse(env4$year==2015,"blue","red")

out<-cbind(groups,comm.dataALLplant,rowSums(comm.dataALLplant[,-1]))
out[order(out$groups),]

#doing my own classification based on Geum and SB dominants
out<-cbind(Sample_name=comm.dataALL$Sample_name,comm.dataALLplant,tot=rowSums(comm.dataALLplant[,-1]))
dm<-out%>%
  filter(GEUROS>0|KOBMYO>0|CARRUP>0,tot>70)
dim(dm)[1]
dmplots<-dm$Sample_name
sb<-out%>%
  filter(GEUROS<1,CARPYR>0|SIBPRO>0|JUNDRU>0)
dim(sb)[1]
sbplots<-sb$Sample_name
dmsbplots<-c(dmplots,sbplots)

comm.dataALL.dm<-comm.dataALL%>%
  filter(Sample_name%in%dmplots)
comm.dataALL.sb<-comm.dataALL%>%
  filter(Sample_name%in%sbplots)

comm.dataALL.dmsb<-rbind(comm.dataALL.dm,comm.dataALL.sb)
lomehidmsb<-rep(c("dm","sb"),each=17)

comm.dataALLplant2<-cbind(lomehi=lomehidmsb,comm.dataALL.dmsb[,4929:4982])
comm.dataALLmicrobe2<-cbind(lomehi=lomehidmsb,comm.dataALL.dmsb[,32:4928])

trts<-as.vector(levels(as.factor(lomehidmsb)))


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
  tempp1<-subset(comm.dataALLplant2, lomehi==trt.temp)
  tempm1<-subset(comm.dataALLmicrobe2, lomehi==trt.temp)
  
  #to make it faster, I should take out zeros here. I am also taking out doubletons and singletons since i will never use them for their correlations and won't use them for correcting qvalue either
  tempp<-cbind(lomehi=tempp1[,1],tempp1[,((which(colSums(tempp1[,2:dim(tempp1)[2]]>0)>2))+1)])
  tempm<-cbind(lomehi=tempm1[,1],tempm1[,((which(colSums(tempm1[,2:dim(tempm1)[2]]>0)>2))+1)])
  
  #in this case the community data started at column 2, so the loop for co-occurrence has to start at that point
  for(b in 2:(dim(tempp)[2])){
    
    results1<-foreach(c=(2):(dim(tempm)[2]),.combine=rbind) %dopar% {
      species1.ab<-sum(tempp[,b])
      species2.ab<-sum(tempm[,c])
      species1.abfreq<-sum(tempp[,b]>0)
      species2.abfreq<-sum(tempm[,c]>0)
      
      if(species1.abfreq>0&species2.abfreq>0){#Ryan had this at 0
        test<-cor.test(tempp[,b],tempm[,c],method="spearman",na.action=na.rm)
        spearmanrho<-test$estimate
        spearmanp.value<-test$p.value
      }    
      if(species1.abfreq <=0 | species2.abfreq <= 0){
        spearmanrho<-0
        spearmanp.value<-1
      }	
      new.row<-c(trts[a],names(tempp)[b],names(tempm)[c],spearmanrho,spearmanp.value,species1.ab,species2.ab,species1.abfreq,species2.abfreq)
    }
    results<-rbind(results,results1)
  }
  print(a/length(trts))
}
head(results)
resultsold<-results
print(Sys.time()-strt) #took 2 minutes
stopCluster(cl)

rownames(resultsold)<-1:dim(resultsold)[1]
results<-data.frame(resultsold[,1:3],(apply(resultsold[,4:dim(resultsold)[2]],2,as.numeric)))
names(results)<-c("trt","taxa1","taxa2","spearmanrho","spearmanp.value","ab1","ab2","ab1freq","ab2freq")
dim(results)



#from run with spearmanm, with a 3 sample frequency cutoff and removal of otus with less than .2% relative abundance
resultsDMSB<-results



###### Edgelist creation #####

#I think the most correct way of doing this is: remove frequencies under whatever limit I want to set (5), do the qvalue correction on all spearman p values (positive and negative rhos), then subset and use only positive rho

edge_listsDMSBa<-cbind(pairs=paste(resultsDMSB$taxa1,resultsDMSB$taxa2),resultsDMSB) #this takes about 10 min
edge_listsDMSBb<-subset(edge_listsDMSBa,ab1freq>5&ab2freq>5) #this takes about 3 min
dim(edge_listsDMSBb)

#take moss out
edge_listsDMSBc<-subset(edge_listsDMSBb,taxa1!="MOSS")

edge_listsDMSBc$qval<-p.adjust(edge_listsDMSBc$spearmanp.value,method="fdr")

#remove negatives
edge_listsDMSB<-subset(edge_listsDMSBc,spearmanrho>0)

#use negatives
edge_listsDMSB<-edge_listsDMSBc
head(edge_listsDMSB)
dim(edge_listsDMSB)
min(edge_listsDMSB$ab1freq)
dim(subset(edge_listsDMSB,trt=="sb"))

unique(subset(edge_listsDMSB,trt=="dm")$taxa1)
unique(subset(edge_listsDMSB,trt=="sb")$taxa1)







#####Plant-microbe interactions#####
#use edge_listsDMSB, cutoff frequency is 4

#DM
inputlo<-subset(edge_listsDMSB,spearmanp.value<.005&abs(spearmanrho)>.65&trt=="dm")[,3:4]
dim(inputlo)
inputlo2<-subset(edge_listsDMSB,spearmanp.value<.005&abs(spearmanrho)>.65&trt=="dm")[,c(3:5)]
graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otu"
colorgraph3<-merge(verticesgraph3,labelsall,"otu",all.y=F,all.x=F,sort=F)
#change all bacteria colors to purple, and plants to green, fungi to pink?
colorgraph3$color[which(colorgraph3$group=="PhotosyntheticBacteria")]<-"#7879BC"
colorgraph3$color[which(colorgraph3$group=="Plant")]<-"#94BA3C"
sizegraph3<-ifelse(colorgraph3$labels=="Plant",8,6)#was 6,4
shapegraph3<-ifelse(colorgraph3$labels=="Plant","square","circle")
linetype3<-ifelse(inputlo2$spearmanrho>0,1,2)
linecolor3<-ifelse(inputlo2$spearmanrho>0,"black","#E95275")
# pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/networkdm.pdf")
plot(graph3,vertex.size=sizegraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=F,vertex.label=NA,vertex.shape=shapegraph3,edge.width=2,edge.color=linecolor3)#,edge.lty=linetype3,,edge.color="gray40"
#dev.off()
sum(linetype3==1) #positive interactions
sum(linetype3==2) #negative interactions
sum(linetype3==2)/(sum(linetype3==1)+sum(linetype3==2))
dim(inputlo)[1]/length(unique(inputlo$taxa1))
dim(inputlo)[1]/17#22

#SB
inputme<-subset(edge_listsDMSB,spearmanp.value<.005&abs(spearmanrho)>.65&trt=="sb")[,3:4]
dim(inputme)
inputme2<-subset(edge_listsDMSB,spearmanp.value<.005&abs(spearmanrho)>.65&trt=="sb")[,3:5]
graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otu"
colorgraph2<-merge(verticesgraph2,labelsall,"otu",all.y=F,all.x=F,sort=F)
#change all bacteria colors to purple, and plants to green, fungi to pink?
colorgraph2$color[which(colorgraph2$group=="PhotosyntheticBacteria")]<-"#7879BC"
colorgraph2$color[which(colorgraph2$group=="Plant")]<-"#94BA3C"
#colorgraph2$color[which(colorgraph2$group=="Fungi")]<-"#E95275"
sizegraph2<-ifelse(colorgraph2$labels=="Plant",8,6)#was 6,4
shapegraph2<-ifelse(colorgraph2$labels=="Plant","square","circle")
linetype2<-ifelse(inputme2$spearmanrho>0,1,2)
linecolor2<-ifelse(inputme2$spearmanrho>0,"black","#E95275")
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Proposals/NSFpreproposal2017/Figs/networksb.pdf")
plot(graph2,vertex.size=sizegraph2,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=F,vertex.label=NA,vertex.shape=shapegraph2,edge.width=2,edge.color=linecolor2)#,edge.color="gray40",edge.lty=linetype2
#dev.off()
sum(linetype2==1) #positive interactions
sum(linetype2==2) #negative interactions
sum(linetype2==2)/(sum(linetype2==1)+sum(linetype2==2))
dim(inputme)[1]/length(unique(inputme$taxa1))
dim(inputme)[1]/7#8

#Sorensons overlap
(2*length(intersect(inputlo$taxa2,inputme$taxa2)))/(length(intersect(inputlo$taxa2,inputme$taxa2))+length(union(inputlo$taxa2,inputme$taxa2)))




