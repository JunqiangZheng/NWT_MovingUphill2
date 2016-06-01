#old code to replace fungi from 18S dataset with fungi from ITS


###### Replace fungi from 18S with ITS ######

datITSr2
datEukr2

head(sample_data(datEukr2))

#figure out which samples are shared across Euk and ITS datasets. Euk doesnt have 81, ITS doesn't have 126 (no reads amplified at all)
temp<-as.data.frame(as.matrix(sample_data(datEukr2)[,"Sample_name"]))
sort(as.numeric(temp$Sample_name))
temp2<-as.data.frame(as.matrix(sample_data(datITSr2)[,"Sample_name"]))
sort(as.numeric(temp2$Sample_name))
cbind(sort(as.numeric(as.character(temp$Sample_name))),sort(as.numeric(as.character(sort(temp2$Sample_name)))))
#temp3<-as.data.frame(as.matrix(sample_data(datBacr2)[,"Sample_name"]))

datEukr3<-subset_samples(datEukr2,Sample_name!="126")
datITSr3<-subset_samples(datITSr2,Sample_name!="81")

#calculated relative abundance of fungi
temp<-otu_table(datEukr3)
temp2<-aggregate.data.frame(temp,by=list(labels=labelsEuk),sum)
rownames(temp2)<-temp2$labels;temp2$labels<-NULL
fungirelabun<-t(temp2["Fungi",])

#remove fungi from euks
datEukr4<-subset_taxa(datEukr3,labels!="Fungi")#labels has no NAs so this works

datEukr4otu<-cbind(sample_data(datEukr4),t(otu_table(datEukr4)))
head(datEukr4otu)[,1:40]
datITSr3otu<-cbind(sample_data(datITSr3),t(otu_table(datITSr3)))
head(datITSr3otu)[,1:40]

#make sample numeric
datEukr4otu$Sample_name<-as.numeric(as.character(datEukr4otu$Sample_name))
datITSr3otu$Sample_name<-as.numeric(as.character(datITSr3otu$Sample_name))

#test if there is denovo name overlap, yes a ton, so relabel the denovos here and in the labels files
nameseuk<-names(datEukr4otu[,-c(1:31)])
namesits<-names(datITSr3otu[,-c(1:31)])
length(nameseuk)
length(namesits)
length(union(nameseuk,namesits))
intersect(nameseuk,namesits)

nameseuk2 <- sub("^", "e", nameseuk )
namesits2 <- sub("^", "i", namesits )

names(datEukr4otu)[-c(1:31)]<-nameseuk2
names(datITSr3otu)[-c(1:31)]<-namesits2

labelsEuk2<-labelsEuk
rownames(labelsEuk2)<-sub("^", "e",rownames(labelsEuk2))
labelsITS2<-labelsITS
rownames(labelsITS2)<-sub("^", "i",rownames(labelsITS2))



#sort things in numeric order by sample
fungirelabun2<-fungirelabun[order(datEukr4otu$Sample_name)]
names(fungirelabun2)<-rownames(fungirelabun)[order(datEukr4otu$Sample_name)]
datEukr5otu<-datEukr4otu[order(datEukr4otu$Sample_name),]
datITSr4otu<-datITSr3otu[order(datITSr3otu$Sample_name),]

#multiply ITS by fungi rel abun
datITSr4otusp<-datITSr4otu[,-c(1:31)]
datITSr4otusp2<-datITSr4otusp*fungirelabun2
rowSums(datITSr4otusp2)#this appears to have worked, the rowSumes (tot rel abun for each sample are the same as the fungirelabun2 file)

#add ITS to the right of the Euk data
datEukr6otu<-cbind(datEukr5otu,datITSr4otusp2)
dim(datEukr6otu)

#check
idenovo44940 in plot 5 is 0.00189401
temp<-otu_table(datITSr2)
head(temp)[,1:35]
temp[which(rownames(temp)=="denovo44940"),]
0.008008008, then scale *0.23651452 (from fungirelabun), checks
