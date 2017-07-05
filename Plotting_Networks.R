#Plotting
#using edge_listsBEP - has fdr correction on all spearman pvalues, and subset to only positive spearman rhos



##### Merge label files #####

#get sequence names from the comm.data files, match with euk or 16S label file, then rbind labelfiles

#Bacteria: photosynthetic phyla are Chloroflexi and Cyanobacteria and some of the Chlorobi (not the Ignavibacteriaceae, there aren't many of these though, the most abundant is present in 7 samples)
labelsBac3<-as.data.frame(labelsBac2)
labelsBac4<-cbind(otu=rownames(labelsBac3),labelsBac3)
labelsBac4$group<-NA

unique(labelsBac4)
head(labelsBac4)
ind<-which(labelsBac4$labels=="Chloroflexi"|labelsBac4$labels=="Cyanobacteria"|labelsBac4$labels=="Chlorobi")
labelsBac4$group[ind]<-"PhotosyntheticBacteria"
labelsBac4$group[-ind]<-"NonphotosyntheticBacteria"

labelsITS2
labelsITS3<-as.data.frame(labelsITS2)
labelsITS4<-cbind(otu=rownames(labelsITS3),labelsITS3)
labelsITS4$group<-"Fungi"

labelsEukN2
labelsEukN3<-as.data.frame(labelsEukN2)
labelsEukN4<-cbind(otu=rownames(labelsEukN3),labelsEukN3)
labelsEukN4$group<-"Metazoa"

labelsEukS2
labelsEukS3<-as.data.frame(labelsEukS2)
labelsEukS4<-cbind(otu=rownames(labelsEukS3),labelsEukS3)
labelsEukS4$group<-NA
ind<-which(labelsEukS4$labels=="Archaeplastida"|labelsEukS4$labels=="Photosynthetic_Stramenopiles"|labelsEukS4$labels=="Photosynthetic_Excavata"|labelsEukS4$labels=="Photosynthetic_Chryptophyceae"|labelsEukS4$labels=="Photosynthetic_Alveolata")
labelsEukS4$group[ind]<-"PhotosyntheticEukaryota"
labelsEukS4$group[-ind]<-"NonphotosyntheticEukaryota"

labelsPlant2<-labelsPlant
labelsPlant2$group<-"Plant"

labels99Met2
labels99Met3<-as.data.frame(labels99Met2)
labels99Met4<-cbind(otu=rownames(labels99Met3),labels99Met3)
labels99Met4$group<-NA
ind<-which(labels99Met4$labels=="Mollusca"|labels99Met4$labels=="Platyhelminthes"|labels99Met4$labels=="Cnidaria"|labels99Met4$labels=="Annelida"|labels99Met4$labels=="Gastrotricha"|labels99Met4$labels=="Nemertea"|labels99Met4$labels=="Entoprocta"|labels99Met4$labels=="Chaetognatha"|labels99Met4$labels=="Brachiopoda")
labels99Met4$group[ind]<-"OtherMetazoa"
ind<-which(labels99Met4$labels=="Nematoda")
labels99Met4$group[ind]<-"Nematoda"
ind<-which(labels99Met4$labels=="Tardigrada")
labels99Met4$group[ind]<-"Tardigrada"
ind<-which(labels99Met4$labels=="Rotifera")
labels99Met4$group[ind]<-"Rotifera"
ind<-which(labels99Met4$labels=="Arthropoda")
labels99Met4$group[ind]<-"Arthopoda"

labels99Nem
labels99Nem2<-cbind(otu=rownames(labels99Nem),labels99Nem)
labels99Nem2$group<-labels99Nem2$labels
labels99Nem2$labels<-"Nematoda"

plantcols<-data.frame(rbind(c("NonphotosyntheticEukaryota","#673482"),
                            c("PhotosyntheticEukaryota","#466D24"),
                            c("Fungi","#F6EC32"),
                            c("Metazoa","#ff9c34"),# ba543d c74e05 #ff9c34 ff9d66
                            c("Plant","#E95275"),#   # 
                            c("PhotosyntheticBacteria","#94BA3C"),
                            c("NonphotosyntheticBacteria","#7879BC"),
                            c("OtherMetazoa","black"),
                            c("Nematoda","gray30"),
                            c("Tardigrada","gray55"),
                            c("Rotifera","gray80"),
                            c("Arthropoda","white"),
                            c("AF","red"),
                            c("AP","red"),
                            c("BF","black"),
                            c("FF","gray30"),
                            c("OM","gray55"),
                            c("PP","gray80"),
                            c("RA","white"),
                            c("unknown","red")))
colnames(plantcols)=c("group","color")

#run this line for all metazoa for 99
#labelsall1<-rbind(labelsBac4,labelsITS4,labelsEukN4,labelsEukS4,labelsPlant2,labels99Met4)

#run this line for only nematodes for 99
labelsall1<-rbind(labelsBac4,labelsITS4,labelsEukN4,labelsEukS4,labelsPlant2,labels99Nem2)

labelsall<-merge(labelsall1,plantcols,"group",all.x=F,all.y=F)
labelsall$color<-as.character(labelsall$color)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/legend.pdf")
plot(c(1,1),c(1,1))
legend("topright",c("Heterotrophic bacteria","Photosynthetic bacteria","Heterotrophic eukaryota","Photosynthetic eukaryota","Metazoa","Fungi","Plants","Bacterial feeder","Fungal feeder","Omnivore","Plant parasite","Root associate"),pt.bg=c("#7879BC","#94BA3C","#673482","#466D24","#ff9c34","#F6EC32","#E95275","black","gray30","gray55","gray80","white"),bty="n",pch=21,cex=1.4)
dev.off()





###### Plotting ######

#All Bac, ITS, Euk small, Euk metazoa 
#Low density
inputlo<-subset(edge_listsBEP,qval<.01&spearmanrho>.6&trt=="lo")[,3:4]
dim(inputlo)
#inputlov<-subset(edge_listsKS32no2b,qval<.05&trt=="lo")
#vertexsizes3<-unique(data.frame(otu=c(as.character(inputlov$taxa1),as.character(inputlov$taxa2)),abun=c(inputlov$ab1,inputlov$ab2)))

graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
#graph2$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph2)
#membership(eb)
#plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otu"
colorgraph3<-merge(verticesgraph3,labelsall,"otu",all.y=F,all.x=F,sort=F)
#sizesgraph3<-merge(verticesgraph3,vertexsizes3,"otu",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuplantbacfuneukf10q.01r.6.pdf")
plot(graph3,vertex.size=4,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2
#dev.off()

#sized by hubs/connectors
sizesgraph3<-ifelse(verticesgraph3$otu%in%hubslo|verticesgraph3$otu%in%connectorslo,8,4)
shapesgraph3<-ifelse(verticesgraph3$otu%in%hubslo,"csquare",'circle')
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuplantbacfuneukf10q.01r.6sizehubcon.pdf")
plot(graph3,vertex.size=sizesgraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA,vertex.shape=shapesgraph3)#vertex.size=log(sizesgraph3$abun)*2   vertex.frame.color=ifelse(verticesgraph3$otu%in%hubslo,"white","black")
dev.off()
#colored by module

colorgraph3b<-membership(cluster_edge_betweenness(graphlo))
colorgraph3c<-data.frame(num=1:54,color=c("#5e91eb",
"#73ba2d",
"#9040b0",
"#3dac31",
"#d452bb",
"#3fd070",
"#b374e7",
"#b2c634",
"#5864d5",
"#d1b524",
"#5567ba",
"#8555a3",
"#6fc65f",
"#c3428c",
"#3b9434",
"#db346f",
"#48ae6a",
"#c93844",
"#64c99c",
"#e65634",
"#43c1c2",
"#b6431f",
"#53a2d5",
"#db7e2d",
"#5569a5",
"#89be4e",
"#d986cc",
"#60912b",
"#aa99df",
"#91a824",
"#d0b847",
"#995584",
"#357429",
"#e36883",
"#3a9674",
"#a0485e",
"#b5bb58",
"#e08eaa",
"#276e4a",
"#e7836b",
"#659258",
"#ab5446",
"#94b86c",
"#9d5723",
"#c7bc7b",
"#52621d",
"#dfa03a",
"#78733a",
"#e1a473",
"#7a8428",
"#ad7b46",
"#aa952d",
"#80621b",
"#a69a55"))
colorgraph3c$color<-as.character(colorgraph3c$color)
colorgraph3d<-data.frame(num=print(colorgraph3b))
colorgraph3d$order<-1:nrow(colorgraph3d)
colorgraph3e<-merge(colorgraph3d,colorgraph3c,"num")
colorgraph3f<-colorgraph3e[order(colorgraph3e$order),]
sizesgraph3<-ifelse(verticesgraph3$otu%in%hubslo|verticesgraph3$otu%in%connectorslo,8,4)
#sizesgraph3<-ifelse(verticesgraph3$otu%in%connectorslo,8,4)
shapesgraph3<-ifelse(verticesgraph3$otu%in%hubslo,"csquare",'circle')
plot(graph3,vertex.size=sizesgraph3,vertex.color=colorgraph3f$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA,vertex.shape=shapesgraph3)#vertex.size=log(sizesgraph3$abun)*2




#Medium density
inputme<-subset(edge_listsBEP,qval<.01&spearmanrho>.6&trt=="me")[,3:4]
dim(inputme)
#inputmev<-subset(edge_listsKS32no2b,qval<.05&trt=="me")
#vertexsizes2<-unique(data.frame(otu=c(as.character(inputmev$taxa1),as.character(inputmev$taxa2)),abun=c(inputmev$ab1,inputmev$ab2)))

graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
#graph2$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph2)
#membership(eb)
#plot(graph2,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otu"
colorgraph2<-merge(verticesgraph2,labelsall,"otu",all.y=F,all.x=F,sort=F)
#sizesgraph2<-merge(verticesgraph2,vertexsizes2,"otuxy",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuplantbacfuneukf10q.01r.6.pdf")
plot(graph2,vertex.size=4,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph2$abun)*2 #dev.off()

#sized by hubs/connectors
sizesgraph2<-ifelse(verticesgraph2$otu%in%hubsme|verticesgraph2$otu%in%connectorsme,8,4)
shapesgraph2<-ifelse(verticesgraph2$otu%in%hubsme,"csquare",'circle')
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuplantbacfuneukf10q.01r.6sizehubcon.pdf")
plot(graph2,vertex.size=sizesgraph2,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA,vertex.shape=shapesgraph2)#vertex.size=log(sizesgraph3$abun)*2   vertex.frame.color=ifelse(verticesgraph3$otu%in%hubslo,"white","black")
dev.off()


#High density
inputhi<-subset(edge_listsBEP,qval<.01&spearmanrho>.6&trt=="hi")[,3:4]
dim(inputhi)
#inputhiv<-subset(edge_listsKS32no2b,qval<.05&trt=="hi")#
#vertexsizes1<-unique(data.frame(otu=c(as.character(inputhiv$taxa1),as.character(inputhiv$taxa2)),abun=c(inputhiv$ab1,inputhiv$ab2)))

graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
#graph1$layout <- layout_in_circle
#eb<-edge.betweenness.community(graph1)
#plot(graph1,vertex.size=2,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otu" #
colorgraph1<-merge(verticesgraph1,labelsall,"otu",all.y=F,all.x=F,sort=F)
#sizesgraph1<-merge(verticesgraph1,vertexsizes1,"otuxy",sort=F)
#sizesgraph1<-ifelse(verticesgraph1$otuxy=="denovo559741",6,2)
sizesgraph1<-ifelse(verticesgraph1$otu%in%hubshi,8,4)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuplantbacfuneukf10q.01r.6.pdf") #f=frequency cutoff 5 included, r=rho cutoff .5
plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2  vertex.label=as.character(colorgraph1$orders)  
#dev.off()

#sized by hubs/connectors
sizesgraph1<-ifelse(verticesgraph1$otu%in%hubshi|verticesgraph1$otu%in%connectorshi,8,4)
shapesgraph1<-ifelse(verticesgraph1$otu%in%hubshi,"csquare",'circle')
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuplantbacfuneukf10q.01r.6sizehubcon.pdf")
plot(graph1,vertex.size=sizesgraph1,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA,vertex.shape=shapesgraph1)#vertex.size=log(sizesgraph3$abun)*2   vertex.frame.color=ifelse(verticesgraph3$otu%in%hubslo,"white","black")
dev.off()


#colored by module
colorgraph1b<-membership(cluster_edge_betweenness(graphhi))
plot(graph1,vertex.size=sizesgraph1,vertex.color=colorgraph1b,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2


eb<-edge.betweenness.community(graph1)
membership(eb)
plot(graph1,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)

eb<-walktrap.community(graph1)
membership(eb)
plot(graph1,vertex.size=4,vertex.color=membership(eb),edge.curved=T,vertex.label=NA)


##if i want to only plot vertices with >1 edge
graph3b<-delete_vertices(graph3, which(degree(graph3)<2))
plot(graph3b,vertex.size=4,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2


temp<-subset(comm.data,greater66plants=="hi")
plot(jitter(temp$Classiculales),temp[,"Unclassified Coccolithales"])
summary(lm(temp$Classiculales~temp[,"Unclassified Coccolithales"]))
cor.test(temp$Classiculales,temp[,"Unclassified Coccolithales"],method="spearman",na.action=na.rm)


taxagraph1<-as.data.frame(get.edgelist(graph1))
colnames(taxagraph1)[1]<-"orders"
taxagraph1a<-unique(merge(taxagraph1,labelfile,"orders",sort=F,all.y=F))
colnames(taxagraph1a)<-c("V1","orders","V1kingdoms")
taxagraph1b<-unique(merge(taxagraph1a,labelfile,"orders",sort=F))
colnames(taxagraph1b)<-c("V2","V1","V1kingdoms","V2kingdoms")
taxagraph1c<-data.frame(V1=taxagraph1b$V1,V2=taxagraph1b$V2,taxagraph1b[,3:4])

taxagraph2<-as.data.frame(get.edgelist(graph2))
colnames(taxagraph2)[1]<-"orders"
taxagraph2a<-unique(merge(taxagraph2,labelfile,"orders",sort=F,all.y=F))
colnames(taxagraph2a)<-c("V1","orders","V1kingdoms")
taxagraph2b<-unique(merge(taxagraph2a,labelfile,"orders",sort=F))
colnames(taxagraph2b)<-c("V2","V1","V1kingdoms","V2kingdoms")
taxagraph2c<-data.frame(V1=taxagraph2b$V1,V2=taxagraph2b$V2,taxagraph2b[,3:4])
















###### Plotting with nematodes ######

#####Nematode-microbe interactions#####
#plot only vertices that are connected to a nematode in all plant densities

#Low density
inputlo<-subset(edge_listsM2g,qval<.05&spearmanrho>.5&trt=="lo")[,3:4]
dim(inputlo)
graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otu"
colorgraph3<-merge(verticesgraph3,labelsall,"otu",all.y=F,all.x=F,sort=F)
sizegraph3<-ifelse(colorgraph3$labels=="Nematoda",8,6)
shapegraph3<-ifelse(colorgraph3$labels=="Nematoda","square","circle")
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuNf4q.05r.5labels.pdf")
plot(graph3,vertex.size=sizegraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.4,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.shape=shapegraph3)# ,vertex.label=NA
title("Low density")
#legend(0,1,c("Heterotrophic bacteria","Photosynthetic bacteria","Heterotrophic eukaryota","Photosynthetic eukaryota","Fungi","Plants","Bacterial feeder","Fungal feeder","Omnivore","Plant parasite","Root associate"),pt.bg=c("#7879BC","#94BA3C","#673482","#466D24","#F6EC32","#E95275","black","gray35","gray55","gray80","white"),bty="n",pch=c(rep(21,6),rep(22,5)),cex=.9)
dev.off()

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/legend2.pdf")
plot(c(1,1),c(1,1))
legend("topright",c("Heterotrophic bacteria","Photosynthetic bacteria","Heterotrophic eukaryota","Photosynthetic eukaryota","Metazoa","Fungi","Plants","Bacterial feeder","Fungal feeder","Omnivore","Plant parasite","Root associate"),pt.bg=c("#7879BC","#94BA3C","#673482","#466D24","#ff9c34","#F6EC32","#E95275","black","gray30","gray55","gray80","white"),bty="n",pch=c(rep(21,7),rep(22,5)),cex=1.4)
dev.off()


#Medium density
inputme<-subset(edge_listsM2g,qval<.05&spearmanrho>.5&trt=="me")[,3:4]
dim(inputme)
graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otu"
colorgraph2<-merge(verticesgraph2,labelsall,"otu",all.y=F,all.x=F,sort=F)
sizegraph2<-ifelse(colorgraph2$labels=="Nematoda",8,6)
shapegraph2<-ifelse(colorgraph2$labels=="Nematoda","square","circle")
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuNf4q.05r.5labels.pdf")
plot(graph2,vertex.size=sizegraph2,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.4,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.shape=shapegraph2)#)#,vertex.label=NA
title("Medium density")
#dev.off()

temp<-subset(comm.dataALLn,lomehi=="me")
cor.test(temp$ndenovo372907,temp$denovo300641,method="spearman",na.action=na.rm)
plot(rank(temp$ndenovo372907),rank(temp$denovo300641))
ndenovo361610


#High density
inputhi<-subset(edge_listsM2g,qval<.05&spearmanrho>.5&trt=="hi")[,3:4]
dim(inputhi)
graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otu" 
colorgraph1<-merge(verticesgraph1,labelsall,"otu",all.y=F,all.x=F,sort=F)
sizegraph1<-ifelse(colorgraph1$labels=="Nematoda",8,6)#was 6,4
shapegraph1<-ifelse(colorgraph1$labels=="Nematoda","square","circle")
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuNf4q.05r.5labels.pdf") #f=frequency cutoff 5 included, r=rho cutoff .5
plot(graph1,vertex.size=sizegraph1,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.3,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.shape=shapegraph1)#,vertex.label=NA
title("High density")
#dev.off()

#network complexity
length(E(graph3))/length(V(graph3))
length(E(graph2))/length(V(graph2))
length(E(graph1))/length(V(graph1))








###### Plotting with Metazoa99 ######

#####Metazoa-microbe interactions#####
#plot only vertices that are connected to a metazoan99 in all plant densities
#use edge_listsM2j

#Low density
inputlo<-subset(edge_listsM2j,qval<.05&spearmanrho>.5&trt=="lo")[,3:4]
dim(inputlo)
graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otu"
colorgraph3<-merge(verticesgraph3,labelsall,"otu",all.y=F,all.x=F,sort=F)
sizegraph3<-ifelse(colorgraph3$labels=="Nematoda",8,6)
shapegraph3<-ifelse(colorgraph3$labels=="Nematoda","square","circle")
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuNf4q.05r.5labels.pdf")
plot(graph3,vertex.size=sizegraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.4,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.shape=shapegraph3)# ,vertex.label=NA
title("Low density")
plot(graph1)
#legend(0,1,c("Heterotrophic bacteria","Photosynthetic bacteria","Heterotrophic eukaryota","Photosynthetic eukaryota","Fungi","Plants","Bacterial feeder","Fungal feeder","Omnivore","Plant parasite","Root associate"),pt.bg=c("#7879BC","#94BA3C","#673482","#466D24","#F6EC32","#E95275","black","gray35","gray55","gray80","white"),bty="n",pch=c(rep(21,6),rep(22,5)),cex=.9)
#dev.off()

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/legend2.pdf")
plot(c(1,1),c(1,1))
legend("topright",c("Heterotrophic bacteria","Photosynthetic bacteria","Heterotrophic eukaryota","Photosynthetic eukaryota","Metazoa","Fungi","Plants","Bacterial feeder","Fungal feeder","Omnivore","Plant parasite","Root associate"),pt.bg=c("#7879BC","#94BA3C","#673482","#466D24","#ff9c34","#F6EC32","#E95275","black","gray30","gray55","gray80","white"),bty="n",pch=c(rep(21,7),rep(22,5)),cex=1.4)
dev.off()


#Medium density
inputme<-subset(edge_listsM2j,qval<.05&spearmanrho>.5&trt=="me")[,3:4]
dim(inputme)
graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otu"
colorgraph2<-merge(verticesgraph2,labelsall,"otu",all.y=F,all.x=F,sort=F)
sizegraph2<-ifelse(colorgraph2$labels=="Nematoda",8,6)
shapegraph2<-ifelse(colorgraph2$labels=="Nematoda","square","circle")
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuNf4q.05r.5labels.pdf")
plot(graph2,vertex.size=sizegraph2,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.4,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.shape=shapegraph2)#)#,vertex.label=NA
title("Medium density")
plot(graph2)
#dev.off()

temp<-subset(comm.dataALLn,lomehi=="me")
cor.test(temp$ndenovo372907,temp$denovo300641,method="spearman",na.action=na.rm)
plot(rank(temp$ndenovo372907),rank(temp$denovo300641))
ndenovo361610


#High density
inputhi<-subset(edge_listsM2j,qval<.05&spearmanrho>.5&trt=="hi")[,3:4]
dim(inputhi)
graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otu" 
colorgraph1<-merge(verticesgraph1,labelsall,"otu",all.y=F,all.x=F,sort=F)
sizegraph1<-ifelse(colorgraph1$labels=="Nematoda",8,6)#was 6,4
shapegraph1<-ifelse(colorgraph1$labels=="Nematoda","square","circle")
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuNf4q.05r.5labels.pdf") #f=frequency cutoff 5 included, r=rho cutoff .5
plot(graph1,vertex.size=sizegraph1,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.3,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.shape=shapegraph1)#,vertex.label=NA
title("High density")
plot(graph1)
#dev.off()







##### Plotting only Bacteria subset ######
#I don't think I ever did this because I can't find anywhere where I created edge_listsBEP16Sc

#Low density
inputlo<-subset(edge_listsBEP16Sc,qval<.01&spearmanrho>.65&trt=="lo")[,3:4]#
dim(inputlo)

graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otuxy"
colorgraph3<-merge(verticesgraph3,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
sizesgraph3<-ifelse(verticesgraph3$otuxy=="denovo559741",5,2)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuplantbaceukf10q.01r.65.pdf")
plot(graph3,vertex.size=4,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2
#dev.off()
#plot(graph3,vertex.size=sizesgraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2


#Medium density
inputme<-subset(edge_listsBEP16Sc,qval<.01&spearmanrho>.65&trt=="me")[,3:4]
dim(inputme)
graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otuxy"
colorgraph2<-merge(verticesgraph2,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#sizesgraph2<-merge(verticesgraph2,vertexsizes2,"otuxy",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuplantbaceukf10q.01r.65.pdf")
plot(graph2,vertex.size=4,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph2$abun)*2
#dev.off()


#High density
inputhi<-subset(edge_listsBEP16Sc,qval<.01&spearmanrho>.65&trt=="hi")[,3:4]
dim(inputhi)
graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otuxy" #change to "order" if doing things by order
colorgraph1<-merge(verticesgraph1,labelsall,"otuxy",all.y=F,all.x=F,sort=F)
#sizesgraph1<-merge(verticesgraph1,vertexsizes1,"otuxy",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuplantbaceukf10q.01r.65.pdf") #f=frequency cutoff 5 included, r=rho cutoff .5
plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2
#dev.off()
#plot(graph1,vertex.size=sizesgraph1,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2









##### Plotting only Eukaryote subset ######
#use edge_listsEuks, make sure the dataset is loaded with the one for the correct frequency cutoff, 5 vs 10

#Low density
inputlo<-subset(edge_listsEuks,qval<.01&spearmanrho>.6&trt=="lo")[,3:4]#
inputlo<-subset(edge_listsEuks,qval<.05&spearmanrho>.5&trt=="lo")[,3:4]#
dim(inputlo)
graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otu"
colorgraph3<-merge(verticesgraph3,labelsall,"otu",all.y=F,all.x=F,sort=F)
#sizesgraph3<-ifelse(verticesgraph3$otuxy=="denovo559741",5,2)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuplantbaceukf10q.01r.65.pdf")
plot(graph3,vertex.size=4,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2
#dev.off()
#plot(graph3,vertex.size=sizesgraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph3$abun)*2


#Medium density
inputme<-subset(edge_listsEuks,qval<.01&spearmanrho>.6&trt=="me")[,3:4]
inputme<-subset(edge_listsEuks,qval<.05&spearmanrho>.5&trt=="me")[,3:4]
dim(inputme)
graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otu"
colorgraph2<-merge(verticesgraph2,labelsall,"otu",all.y=F,all.x=F,sort=F)
#sizesgraph2<-merge(verticesgraph2,vertexsizes2,"otuxy",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuplantbaceukf10q.01r.65.pdf")
plot(graph2,vertex.size=4,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#vertex.size=log(sizesgraph2$abun)*2
#dev.off()


#High density
inputhi<-subset(edge_listsEuks,qval<.01&spearmanrho>.6&trt=="hi")[,3:4]
inputhi<-subset(edge_listsEuks,qval<.05&spearmanrho>.5&trt=="hi")[,3:4]
dim(inputhi)
graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otu" #change to "order" if doing things by order
colorgraph1<-merge(verticesgraph1,labelsall,"otu",all.y=F,all.x=F,sort=F)
#sizesgraph1<-merge(verticesgraph1,vertexsizes1,"otuxy",sort=F)
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuplantbaceukf10q.01r.65.pdf") #f=frequency cutoff 5 included, r=rho cutoff .5
plot(graph1,vertex.size=4,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2
#dev.off()
#plot(graph1,vertex.size=sizesgraph1,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA)#,vertex.size=log(sizesgraph1$abun)*2


#network complexity
length(E(graph3))/length(V(graph3))
length(E(graph2))/length(V(graph2))
length(E(graph1))/length(V(graph1))




#####Plant-microbe interactions#####
#use edge_listsPlants, cutoff frequency is 5
#extract only vertices that are connected to a plant in all plant densities to determine if plant correlations are still there
#figure out why there is more richness in high plant density but the same number of network vertices - is the added diversity all low abundance?

#Low density
inputlo<-subset(edge_listsPlants,qval<.05&spearmanrho>.5&trt=="lo")[,3:4]
dim(inputlo)
#inputlo2<-inputlo[which(inputlo$taxa1%in%plantlabels$otu|inputlo$taxa2%in%plantlabels$otu),]
graph3<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
verticesgraph3<-as.data.frame(rownames(as.matrix(V(graph3))))
colnames(verticesgraph3)<-"otu"
colorgraph3<-merge(verticesgraph3,labelsall,"otu",all.y=F,all.x=F,sort=F)
sizegraph3<-ifelse(colorgraph3$labels=="Plant",8,6)#was 6,4
shapegraph3<-ifelse(colorgraph3$labels=="Plant","square","circle")
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/lodensityotuplantf5q.05r.5.pdf")
plot(graph3,vertex.size=sizegraph3,vertex.color=colorgraph3$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA,vertex.shape=shapegraph3)#
#dev.off()


#Medium density
inputme<-subset(edge_listsPlants,qval<.05&spearmanrho>.5&trt=="me")[,3:4]
dim(inputme)
graph2<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
verticesgraph2<-as.data.frame(rownames(as.matrix(V(graph2))))
colnames(verticesgraph2)<-"otu"
colorgraph2<-merge(verticesgraph2,labelsall,"otu",all.y=F,all.x=F,sort=F)
sizegraph2<-ifelse(colorgraph2$labels=="Plant",8,6)#was 6,4
shapegraph2<-ifelse(colorgraph2$labels=="Plant","square","circle")
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/medensityotuplantf5q.05r.5.pdf")
plot(graph2,vertex.size=sizegraph2,vertex.color=colorgraph2$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.label=NA,vertex.shape=shapegraph2)
#dev.off()

#High density
inputhi<-subset(edge_listsPlants,qval<.05&spearmanrho>.5&trt=="hi")[,3:4]
dim(inputhi)
graph1<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))
verticesgraph1<-as.data.frame(rownames(as.matrix(V(graph1))))
colnames(verticesgraph1)<-"otu"
colorgraph1<-merge(verticesgraph1,labelsall,"otu",all.y=F,all.x=F,sort=F)
sizegraph1<-ifelse(colorgraph1$labels=="Plant",8,6)#was 6,4
shapegraph1<-ifelse(colorgraph1$labels=="Plant","square","circle")
#pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/hidensityotuplantf5q.05r.5.pdf") #f=frequency cutoff 5 included, r=rho cutoff .5
plot(graph1,vertex.size=sizegraph1,vertex.color=colorgraph1$color,vertex.label.cex=.8,vertex.label.dist=.1,vertex.label.color="black",edge.curved=T,edge.color="gray40",vertex.shape=shapegraph1,vertex.label=NA)#,vertex.label=NA
#dev.off()








