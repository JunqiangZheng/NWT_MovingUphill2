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
length(E(graphlo))/length(V(graphlo))
length(E(graphme))/length(V(graphme))
length(E(graphhi))/length(V(graphhi))
graph.density(graphlo)
graph.density(graphme)
graph.density(graphhi)
transitivity(graphlo, type="global")
transitivity(graphme, type="global")
transitivity(graphhi, type="global")#clustering, triangle
modularity(graphlo, membership(walktrap.community(graphlo)))
modularity(graphme, membership(walktrap.community(graphme)))
modularity(graphhi, membership(walktrap.community(graphhi)))

#link denstiy, average number of links, not necessary
graphlo.adj<-get.adjacency(graphlo,sparse=F,type="upper")# in older versions of igraph the default was sparse=F, but now you must specify, other wise you get a matrix of 1s and .s. if you don't use type=upper some metrics are multiplied by 2 because it is using the whole matrix which has double the info (not just upper or lower)
graphlo.adj.properties<-GenInd(graphlo.adj)
graphlo.adj.properties$N
graphlo.adj.properties$Ltot
graphlo.adj.properties$LD   

#I should look into the clustering method more
modularity(graphlo, membership(cluster_edge_betweenness(graphlo)))
modularity(graphme, membership(cluster_edge_betweenness(graphme)))
modularity(graphhi, membership(cluster_edge_betweenness(graphhi)))

#checking graph density values, this checks, not sure why the person in stack overflow had problems with it.
#density = mean degree / (n-1)
mean(degree(graphme))/(vcount(graphme)-1)

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

statslo2<-statslo[order(statslo$norm_degree,decreasing=T),]
hist(statslo2$norm_degree) #highly non-normal, can't use standard deviation to get at "outliers"
statshi2<-statshi[order(statshi$norm_degree,decreasing=T),]
hist(statshi2$norm_degree)


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



#Calculating Zi and Pi by "hand"
#info from bipartite pdf
#any vertex linked with only one other vertex in a module will have an NA for zi b/c the sd=0 (or any vertex linked with the same number of taxa in each module)
p = 1 - sum( (k.it/k.i)^2) # among-module connectivity = participation coefficient P in Guimerà & Amaral
z = (k.is - ks.bar) / SD.ks # within-module degree
k.is = number of links of i to other species in its own module s
ks.bar = average k.is of all species in module  
SD.ks = standard deviation of k.is of all species in module s
k.it = number of links of species i to module t
k.i = degree of species i
#Note that for any species alone (in its level) in a module the z-value will be NaN, since then SD.ks is 0. This is a limitation of the way the z-value is defined (in multiples of degree/strength standard deviations).
#Olesen et al. (2006) give critical c and z values of 0.62 and 2.5, respectively. Species exceeding these values are deemed connectors or hubs of a network. The justification of these thresholds remains unclear to me. They may also not apply for the quantitative version.

#assuming the same modules as calculated above in "modularity"
graphlo #has vertex links

#graphlo Zi
graphlomem<-membership(cluster_edge_betweenness(graphlo)) #621 vertices
graphloedges<-cbind(as.character(inputlov[,3]),as.character(inputlov[,4]))

zipioutputlo<-data.frame(otu=rownames(as.matrix(V(graphlo))),zi=rep(NA,length(graphlomem)))

for(i in 1:max(graphlomem)){
  #extract the vertices for each module
  graphlovertices<-names(graphlomem)[which(graphlomem==i)]
  
  #subset only those vertices from the network
  #E(graphlo)
  #inputlov[,3:4] the dataframe with the two pairs as columns
  ind<-which(graphloedges[,1]%in%graphlovertices&graphloedges[,2]%in%graphlovertices)
  newgraphinput<-graphloedges[ind,]
  #newgraph<-simplify(graph.edgelist(as.matrix(newgraphinput),directed=FALSE))
  newgraph<-simplify(graph.edgelist(matrix(newgraphinput,ncol=2),directed=FALSE))
  
  #calculate zi for each vertex in the module
  k.is<-degree(newgraph,normalized=F)
  ks.bar<-mean(k.is)
  SD.ks<-sd(k.is)
  z <- (k.is - ks.bar) / SD.ks
  
  #put in correct slot in output file
  ind<-match(names(z),zipioutputlo$otu)#ind <- ind[!is.na(ind)]
  zipioutputlo$zi[ind]<-z
}

sort( zipioutputlo$zi)
hist( zipioutputlo$zi)


#graphlo pi

zipioutputlo$pi<-NA
zipioutputlo$degree<-NA


for(i in 1:length(V(graphlo))){
  #get a node
  nodei<-names(graphlomem[i])
  memi<-graphlomem[i]

  #how many/which different nodes is it connected to
  ind<-which(graphloedges[,1]%in%nodei|graphloedges[,2]%in%nodei)#matches degree(graphlo)[i]
  connectorsi<-c(graphloedges[ind,1],graphloedges[ind,2])
  connectorsi<-connectorsi[connectorsi!=nodei]
  
  #how many different modules is it connected to
  #extract modules for each connector
  temp<-print(graphlomem) #need to make this not actually print capture.output does weird things with rownmaes
  ind<-match(connectorsi,names(temp))
  connectorsi.m<-data.frame(graphlomem[ind])
  
  #extract only connectors of a different module
  #connectorsi.m2<-data.frame(connectorsi.m[which(connectorsi.m!=memi)])
  connectorsi.m$ones<-1
  colnames(connectorsi.m)[1]<-"membership"
  k.it<-aggregate.data.frame((connectorsi.m$ones),by=list(membership=connectorsi.m$membership),sum)$x
  
  #calculate ci
  k.i<-degree(graphlo)[i]
  p.i = 1 - sum( (k.it/k.i)^2)
  zipioutputlo$pi[i]<-p.i
  zipioutputlo$degree[i]<-k.i
}

sort(zipioutputlo$pi)
hist(zipioutputlo$pi)
zipioutputlo[order(zipioutputlo$pi,decreasing=T),]




#Functionalizing zi and pi calculations
#assuming the same modules as calculated above in "modularity"

zipi<-function(inputfile,graphfile){

  # Zi
  graphlomem<-membership(cluster_edge_betweenness(graphfile)) #621 vertices
  graphloedges<-cbind(as.character(inputfile[,3]),as.character(inputfile[,4]))
  
  zipioutput<-data.frame(otu=rownames(as.matrix(V(graphfile))),zi=rep(NA,length(graphlomem)))
  
  for(i in 1:max(graphlomem)){
    #extract the vertices for each module
    graphlovertices<-names(graphlomem)[which(graphlomem==i)]
    
    #subset only those vertices from the network
    ind<-which(graphloedges[,1]%in%graphlovertices&graphloedges[,2]%in%graphlovertices)
    newgraphinput<-graphloedges[ind,]
    newgraph<-simplify(graph.edgelist(matrix(newgraphinput,ncol=2),directed=FALSE))
    
    #calculate zi for each vertex in the module
    k.is<-degree(newgraph,normalized=F)
    ks.bar<-mean(k.is)
    SD.ks<-sd(k.is)
    z <- (k.is - ks.bar) / SD.ks
    
    #put in correct slot in output file
    ind<-match(names(z),zipioutput$otu)#ind <- ind[!is.na(ind)]
    zipioutput$zi[ind]<-z
  }

  #Pi
  zipioutput$pi<-NA
  zipioutput$degree<-NA

  for(i in 1:length(V(graphfile))){
    #get a node
    nodei<-names(graphlomem[i])
    memi<-graphlomem[i]
    
    #how many/which different nodes is it connected to
    ind<-which(graphloedges[,1]%in%nodei|graphloedges[,2]%in%nodei)#matches degree(graphlo)[i]
    connectorsi<-c(graphloedges[ind,1],graphloedges[ind,2])
    connectorsi<-connectorsi[connectorsi!=nodei]
    
    #how many different modules is it connected to
    #extract modules for each connector
    temp<-print(graphlomem) #need to make this not actually print capture.output does weird things with rownmaes
    ind<-match(connectorsi,names(temp))
    connectorsi.m<-data.frame(graphlomem[ind])
    
    #extract only connectors of a different module
    #connectorsi.m2<-data.frame(connectorsi.m[which(connectorsi.m!=memi)])
    connectorsi.m$ones<-1
    colnames(connectorsi.m)[1]<-"membership"
    k.it<-aggregate.data.frame((connectorsi.m$ones),by=list(membership=connectorsi.m$membership),sum)$x
    
    #calculate ci
    k.i<-degree(graphfile)[i]
    p.i = 1 - sum( (k.it/k.i)^2)
    zipioutput$pi[i]<-p.i
    zipioutput$degree[i]<-k.i
  }
  return(zipioutput)
}

zipioutputlo<-zipi(inputlov,graphlo)
zipioutputme<-zipi(inputmev,graphme)
zipioutputhi<-zipi(inputhiv,graphhi)

plot(zipioutputlo$pi,zipioutputlo$zi)
abline(h=2.5)
abline(v=.62)
plot(zipioutputme$pi,zipioutputme$zi)
abline(h=2.5)
abline(v=.62)
plot(zipioutputhi$pi,zipioutputhi$zi)
abline(h=2.5)
abline(v=.62)

hubslo<-zipioutputlo$otu[which(zipioutputlo$zi>2.5)]
hubsme<-zipioutputme$otu[which(zipioutputme$zi>2.5)]
hubshi<-zipioutputhi$otu[which(zipioutputhi$zi>2.5)]

which(zipioutputlo$zi>2.5)
which(zipioutputme$zi>2.5)
which(zipioutputhi$zi>2.5)
which(zipioutputlo$pi>.62)
which(zipioutputme$pi>.62)
which(zipioutputhi$pi>.62)

zipioutputlo[order(zipioutputlo$degree,decreasing=T),]

zipioutputhi[which(zipioutputhi$zi>2.5),]





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






#Intersection of edges of two networks
graph.intersection(graphlo,graphme)#union of the vertex names and an intersection of the edges
graph.union(graphlo,graphme)
#jaccard type similarity: intersection/union of edges
96/3791
#sorensons type similarity: 2*intersection/(sum of number of edges in each), I think I will go with sorensons it is more intuitive because it gives a value of .5 if half of the interactions in group 1 are shared (whereas jaccard gives a value of .25)
2*96/(2077+1810)=0.049395

graph.intersection(graphlo,graphhi)
graph.union(graphlo,graphhi)
2*43/(2077+1441)=0.024446

graph.intersection(graphme,graphhi)
graph.union(graphme,graphhi)
2*47/(1810+1441)=0.02891

#intersection across all three plant density levels
temp<-graph.intersection(graphlo,graphme)
graph.intersection(temp,graphhi)
3*9/(2077+1810+1441)=0.0050676

#for nematode networks
inputlo<-subset(edge_listsM2g,qval<.05&spearmanrho>.5&trt=="lo")[,3:4]
graphlon<-simplify(graph.edgelist(as.matrix(inputlo),directed=FALSE))
inputme<-subset(edge_listsM2g,qval<.05&spearmanrho>.5&trt=="me")[,3:4]
graphmen<-simplify(graph.edgelist(as.matrix(inputme),directed=FALSE))
inputhi<-subset(edge_listsM2g,qval<.05&spearmanrho>.5&trt=="hi")[,3:4]
graphhin<-simplify(graph.edgelist(as.matrix(inputhi),directed=FALSE))

graph.intersection(graphlon,graphmen)
graph.intersection(graphmen,graphhin)
graph.intersection(graphlon,graphhin)
#graph.union(graphlon,graphmen)



#Intersection of vertices
length(intersect(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphlo))))) #283
length(union(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphlo)))))
V(graphlo)
V(graphme)
2*283/(615+621)=0.4579

length(intersect(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphhi)))))
length(union(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphhi)))))
V(graphme)
V(graphhi)
2*275/(615+716)=0.41322

length(intersect(row.names(as.matrix(V(graphhi))),row.names(as.matrix(V(graphlo)))))
length(union(row.names(as.matrix(V(graphhi))),row.names(as.matrix(V(graphlo)))))
V(graphlo)
V(graphhi)
2*228/(621+716)=0.34106

temp<-intersect(row.names(as.matrix(V(graphme))),row.names(as.matrix(V(graphlo))))
intersect(temp,row.names(as.matrix(V(graphhi))))
3*153/(615+621+716)=0.2351434

#Taxa in low medium high networks
row.names(as.matrix(V(graphlo)))
row.names(as.matrix(V(graphme))) #SENFRE, MOSS
row.names(as.matrix(V(graphhi))) #ANGGRA FESBRA CIRSCO


#Nematode networks
length(intersect(row.names(as.matrix(V(graphmen))),row.names(as.matrix(V(graphlon))))) #1 nematode
intersect(row.names(as.matrix(V(graphmen))),row.names(as.matrix(V(graphlon))))
V(graphlon)
V(graphmen)
2*1/(12+23)=0.05714286

mdenovo281318
trophicgroup[trophicgroup$X.OTU.ID=="denovo281318",]

length(intersect(row.names(as.matrix(V(graphmen))),row.names(as.matrix(V(graphhin))))) #4 (3 nematodes and 1 bacteria)
intersect(row.names(as.matrix(V(graphmen))),row.names(as.matrix(V(graphhin))))
V(graphmen)
V(graphhin)
2*4/(23+86)=0.0733945

length(intersect(row.names(as.matrix(V(graphlon))),row.names(as.matrix(V(graphhin))))) #3 nematodes
intersect(row.names(as.matrix(V(graphlon))),row.names(as.matrix(V(graphhin))))
V(graphlon)
V(graphhin)
2*3/(12+86)=0.061224

temp<-intersect(row.names(as.matrix(V(graphmen))),row.names(as.matrix(V(graphlon))))
intersect(temp,row.names(as.matrix(V(graphhin))))
3*1/(12+23+86)


#Taxa in low me hi nematode networks
temp<-row.names(as.matrix(V(graphlon)))
trophicgroup[trophicgroup$otu2%in%temp,]

temp<-row.names(as.matrix(V(graphmen)))
trophicgroup[trophicgroup$otu2%in%temp,]

temp<-row.names(as.matrix(V(graphhin)))
trophicgroup[trophicgroup$otu2%in%temp,]

#how many microbes they interacted wtih
temp<-data.frame(otu2=sort(c(as.character(inputlo$taxa1),as.character(inputlo$taxa2))))
temp2<-count(temp,otu2) #in plyr package
merge(temp2,trophicgroup,"otu2")

temp<-data.frame(otu2=sort(c(as.character(inputme$taxa1),as.character(inputme$taxa2))))
temp2<-count(temp,otu2)
merge(temp2,trophicgroup,"otu2")

temp<-data.frame(otu2=sort(c(as.character(inputhi$taxa1),as.character(inputhi$taxa2))))
temp2<-count(temp,otu2)
merge(temp2,trophicgroup,"otu2")
