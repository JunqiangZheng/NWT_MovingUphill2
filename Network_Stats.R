#NetworkStats

#looking at p value and rho vs frequency
inputlov<-subset(edge_listsBEP,qval<.02&spearmanrho>.6)
min(inputlov$spearmanrho)
inputlov2<-subset(edge_listsBEP,qval<.01&spearmanrho>.7)
mean(c(inputlov$ab1freq,inputlov$ab2freq))
mean(c(inputlov2$ab1freq,inputlov2$ab2freq))
temp<-((inputlov$ab1freq+inputlov$ab2freq)/2)
plot((inputlov$ab1freq+inputlov$ab2freq)/2,inputlov$spearmanrho)
abline(lm(inputlov$spearmanrho~temp))
summary(lm(inputlov$spearmanrho~temp))
plot(inputlov$qval,inputlov$spearmanrho)
#pairs with a higher freq, have higher pvalues and lower rhos: making the cutoff more stringent you bias towards infrequent taxa

inputlov<-subset(edge_listsBEP,qval<.01&spearmanrho>.6&trt=="lo")
inputmev<-subset(edge_listsBEP,qval<.01&spearmanrho>.6&trt=="me")
inputhiv<-subset(edge_listsBEP,qval<.01&spearmanrho>.6&trt=="hi")
#dim(inputlov)
#dim(inputmev)
#dim(inputhiv)

graphlo<-simplify(graph.edgelist(as.matrix(inputlov[,3:4]),directed=FALSE))
graphme<-simplify(graph.edgelist(as.matrix(inputmev[,3:4]),directed=FALSE))
graphhi<-simplify(graph.edgelist(as.matrix(inputhiv[,3:4]),directed=FALSE))

length(V(graphlo))#
length(V(graphme))#
length(V(graphhi))#
length(E(graphlo))#
length(E(graphme))#
length(E(graphhi))#
graph.density(graphlo)
graph.density(graphme)
graph.density(graphhi)
transitivity(graphlo, type="global")
transitivity(graphme, type="global")
transitivity(graphhi, type="global")#clustering, triangle
modularity(graphlo, membership(walktrap.community(graphlo)))
modularity(graphme, membership(walktrap.community(graphme)))
modularity(graphhi, membership(walktrap.community(graphhi)))

#I should look into the clustering method more
modularity(graphlo, membership(cluster_edge_betweenness(graphlo)))
modularity(graphme, membership(cluster_edge_betweenness(graphme)))
modularity(graphhi, membership(cluster_edge_betweenness(graphhi)))




#get random transitivity values for hi me lo, this is just one random draw, could do many
transitivity(erdos.renyi.game(length(V(graphlo)),length(E(graphlo)),type="gnm"))
transitivity(erdos.renyi.game(length(V(graphme)),length(E(graphme)),type="gnm"))
transitivity(erdos.renyi.game(length(V(graphhi)),length(E(graphhi)),type="gnm"))

statslo<-data.frame(otu=row.names((as.matrix(degree(graphlo,normalized=TRUE)))),norm_degree=(as.matrix(degree(graphlo,normalized=TRUE))),closeness=(as.matrix(closeness(graphlo,normalized=TRUE))),betweenness=(as.matrix(betweenness(graphlo,normalized=TRUE))))
statsme<-data.frame(otu=row.names((as.matrix(degree(graphme,normalized=TRUE)))),norm_degree=(as.matrix(degree(graphme,normalized=TRUE))),closeness=(as.matrix(closeness(graphme,normalized=TRUE))),betweenness=(as.matrix(betweenness(graphme,normalized=TRUE))))
statshi<-data.frame(otu=row.names((as.matrix(degree(graphhi,normalized=TRUE)))),norm_degree=(as.matrix(degree(graphhi,normalized=TRUE))),closeness=(as.matrix(closeness(graphhi,normalized=TRUE))),betweenness=(as.matrix(betweenness(graphhi,normalized=TRUE))))


#Normalized degree
statslo1<-statslo[order(statslo$norm_degree,decreasing=T),][1:5,]
statsme1<-statsme[order(statsme$norm_degree,decreasing=T),][1:5,]
statshi1<-statshi[order(statshi$norm_degree,decreasing=T),][1:5,]

tax_table(datBac3)[which(rownames(tax_table(datBac3))%in%substring(statslo1$otu,2)),]
tax_table(datBac3)[which(rownames(tax_table(datBac3))%in%substring(statsme1$otu,2)),]
tax_table(datBac3)[which(rownames(tax_table(datBac3))%in%substring(statshi1$otu,2)),]


#Betweenness centrality
statslo1<-statslo[order(statslo$betweenness,decreasing=T),][1:5,]
statsme1<-statsme[order(statsme$betweenness,decreasing=T),][1:5,]
statshi1<-statshi[order(statshi$betweenness,decreasing=T),][1:5,]

labelsall[labelsall$otu%in%rownames(statslo1),]
#all are bacteria
tax_table(datBac3)[which(rownames(tax_table(datBac3))%in%substring(statslo1$otu,2)),]

labelsall[labelsall$otu%in%rownames(statsme1),]
#all are bacteria
tax_table(otufile16S)[which(rownames(tax_table(otufile16S))%in%statsme1$otu),]

labelsall[labelsall$otu%in%rownames(statshi1),]
#all are bacteria
tax_table(otufile16S)[which(rownames(tax_table(otufile16S))%in%statshi1$otu),]

labelsall[labelsall$otu=="denovo143776",]
tax_table(dats2)[rownames(tax_table(dats2))=="denovo143776",]#euk
tax_table(dat16Ss2)[rownames(tax_table(dat16Ss2))=="denovo462276",]



# Number of plants bacteria and euks in networks
verticesgraphlo<-as.data.frame(rownames(as.matrix(V(graphlo))))
verticesgraphlo2<-labelsall[which(labelsall$otu%in%verticesgraphlo[,1]),]
length(which(verticesgraphlo2$group=="NonphotosyntheticBacteria"|verticesgraphlo2$group=="PhotosyntheticBacteria"))
length(which(verticesgraphlo2$group=="Fungi"))
length(which(verticesgraphlo2$group=="PhotosyntheticEukaryota"|verticesgraphlo2$group=="NonphotosyntheticEukaryota"))
length(which(verticesgraphlo2$group=="Plant"))

verticesgraphme<-as.data.frame(rownames(as.matrix(V(graphme))))
verticesgraphme2<-labelsall[which(labelsall$otu%in%verticesgraphme[,1]),]
length(which(verticesgraphme2$group=="NonphotosyntheticBacteria"|verticesgraphme2$group=="PhotosyntheticBacteria"))
length(which(verticesgraphme2$group=="Fungi"))
length(which(verticesgraphme2$group=="PhotosyntheticEukaryota"|verticesgraphme2$group=="NonphotosyntheticEukaryota"))
length(which(verticesgraphme2$group=="Plant"))

verticesgraphhi<-as.data.frame(rownames(as.matrix(V(graphhi))))
verticesgraphhi2<-labelsall[which(labelsall$otu%in%verticesgraphhi[,1]),]
length(which(verticesgraphhi2$group=="NonphotosyntheticBacteria"|verticesgraphhi2$group=="PhotosyntheticBacteria"))
length(which(verticesgraphhi2$group=="Fungi"))
length(which(verticesgraphhi2$group=="PhotosyntheticEukaryota"|verticesgraphhi2$group=="NonphotosyntheticEukaryota"))
length(which(verticesgraphhi2$group=="Plant"))
length(which(verticesgraphhi2$group=="Metazoa"))

#intersection of edges of two graphs
graph.intersection(graphlo,graphme)#union of the vertex names and an intersection of the edges
graph.union(graphlo,graphme)

graph.intersection(graphlo,graphhi)
graph.union(graphlo,graphhi)

graph.intersection(graphme,graphhi)
graph.union(graphme,graphhi)

#for plant only network (run the graph in the Plotting_Networks scripts)
graph.intersection(graph1,graph2)
graph.union(graph1,graph2)


#intersection of vertices
length(intersect(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphlo)))))
length(union(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphlo)))))

length(intersect(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphhi)))))
length(union(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphhi)))))

length(intersect(row.names(as.matrix(V(graphhi))),row.names(as.matrix(V(graphlo)))))
length(union(row.names(as.matrix(V(graphhi))),row.names(as.matrix(V(graphlo)))))


tax_table(dats2)[rownames(tax_table(dats2))=="denovo215842",]
tax_table(dats2)[rownames(tax_table(dats2))=="denovo12769",]






