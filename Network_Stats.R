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




######  Functionalizing zi and pi calculations ######
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
zipioutputlo$density<-"Low"
zipioutputme<-zipi(inputmev,graphme)
zipioutputme$density<-"Medium"
zipioutputhi<-zipi(inputhiv,graphhi)
zipioutputhi$density<-"High"

plot(zipioutputlo$pi,zipioutputlo$zi)
abline(h=2.5)
abline(v=.62)
points(zipioutputme$pi,zipioutputme$zi, col=2)
abline(h=2.5)
abline(v=.62)
points(zipioutputhi$pi,zipioutputhi$zi,col=3)
abline(h=2.5)
abline(v=.62)

zipioutputall<-rbind(zipioutputlo,zipioutputme,zipioutputhi)

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/zipi",width=7, height=4.7) #for two panels width=7, height=3.5)
ggplot(zipioutputall,aes(x=pi,y=zi,color=density))+# as.numeric(fert),color=species
  labs(x="Among-module connectivity (Pi)",y="Within-module connectivity (Zi)")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  xlim(0,1)+
  ylim(-2,6)+
  geom_point(size=1.4)+
  geom_hline(yintercept=2.5)+
  geom_vline(xintercept=.62)
  #geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  #facet_wrap(~type,scales="free")
dev.off()


hubslo<-zipioutputlo$otu[which(zipioutputlo$zi>2.5)]
hubsme<-zipioutputme$otu[which(zipioutputme$zi>2.5)]
hubshi<-zipioutputhi$otu[which(zipioutputhi$zi>2.5)]

connectorslo<-zipioutputlo$otu[which(zipioutputlo$pi>.62)]
connectorsme<-zipioutputme$otu[which(zipioutputme$pi>.62)]
connectorshi<-zipioutputhi$otu[which(zipioutputhi$pi>.62)]

which(zipioutputlo$zi>2.5)
which(zipioutputme$zi>2.5)
which(zipioutputhi$zi>2.5)
which(zipioutputlo$pi>.62)
which(zipioutputme$pi>.62)
which(zipioutputhi$pi>.62)

zipioutputlo[order(zipioutputlo$degree,decreasing=T),]

#Zi: extracting names of taxa for lo me hi

#Low
zlobac<-zipioutputlo[which(zipioutputlo$zi>2.5),]
zlobac2<-zlobac[which(substring(zlobac$otu,1,1)=="b"),]
zlobac2$otu2<-substring(zlobac2$otu,2)

zlofun<-zipioutputlo[which(zipioutputlo$zi>2.5),]
zlofun2<-zlofun[which(substring(zlofun$otu,1,1)=="i"),]
zlofun2$otu2<-substring(zlofun2$otu,2)

zlobac3<-data.frame(tax_table(datBac3)[which(rownames(tax_table(datBac3))%in%substring(zlobac2$otu,2)),])
zlobac3$otu2<-rownames(zlobac3)
zlobac4<-merge(zlobac3,zlobac2,"otu2")

zlofun3<-data.frame(tax_table(datITS3)[which(rownames(tax_table(datITS3))%in%substring(zlofun2$otu,2)),])
zlofun3$otu2<-rownames(zlofun3)
zlofun4<-merge(zlofun3,zlofun2,"otu2")

write.csv(zlobac4,"~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/zlobac4.csv")

#Med, only bacteria
zmebac<-zipioutputme[which(zipioutputme$zi>2.5),]
zmebac2<-zmebac[which(substring(zmebac$otu,1,1)=="b"),]
zmebac2$otu2<-substring(zmebac2$otu,2)

zmebac3<-data.frame(tax_table(datBac3)[which(rownames(tax_table(datBac3))%in%substring(zmebac2$otu,2)),])
zmebac3$otu2<-rownames(zmebac3)
zmebac4<-merge(zmebac3,zmebac2,"otu2")

zmebac4[order(zmebac4$zi,decreasing=T),]

#Hi, only bacteria
zhibac<-zipioutputhi[which(zipioutputhi$zi>2.5),]
zhibac2<-zhibac[which(substring(zhibac$otu,1,1)=="b"),]
zhibac2$otu2<-substring(zhibac2$otu,2)

zhibac3<-data.frame(tax_table(datBac3)[which(rownames(tax_table(datBac3))%in%substring(zhibac2$otu,2)),])
zhibac3$otu2<-rownames(zhibac3)
zhibac4<-merge(zhibac3,zhibac2,"otu2")

zhibac4[order(zhibac4$zi,decreasing=T),]


#Pi:

#Low, only bacteria
plobac<-zipioutputlo[which(zipioutputlo$pi>.62),]
plobac2<-plobac[which(substring(plobac$otu,1,1)=="b"),]
plobac2$otu2<-substring(plobac2$otu,2)

plobac3<-data.frame(tax_table(datBac3)[which(rownames(tax_table(datBac3))%in%substring(plobac2$otu,2)),])
plobac3$otu2<-rownames(plobac3)
plobac4<-merge(plobac3,plobac2,"otu2")

plobac4[order(plobac4$pi,decreasing=T),]

#Med, bacteria and fungi
pmebac<-zipioutputme[which(zipioutputme$pi>.62),]
pmebac2<-pmebac[which(substring(pmebac$otu,1,1)=="b"),]
pmebac2$otu2<-substring(pmebac2$otu,2)

pmebac3<-data.frame(tax_table(datBac3)[which(rownames(tax_table(datBac3))%in%substring(pmebac2$otu,2)),])
pmebac3$otu2<-rownames(pmebac3)
pmebac4<-merge(pmebac3,pmebac2,"otu2")

pmefun<-zipioutputme[which(zipioutputme$pi>.62),]
pmefun2<-pmefun[which(substring(pmefun$otu,1,1)=="i"),]
pmefun2$otu2<-substring(pmefun2$otu,2)

pmefun3<-data.frame(tax_table(datITS3)[which(rownames(tax_table(datITS3))%in%substring(pmefun2$otu,2)),])
pmefun3$otu2<-rownames(pmefun3)
pmefun4<-merge(pmefun3,pmefun2,"otu2")


#Hi, only bacteria
phibac<-zipioutputhi[which(zipioutputhi$pi>.62),]
phibac2<-phibac[which(substring(phibac$otu,1,1)=="b"),]
phibac2$otu2<-substring(phibac2$otu,2)

phibac3<-data.frame(tax_table(datBac3)[which(rownames(tax_table(datBac3))%in%substring(phibac2$otu,2)),])
phibac3$otu2<-rownames(phibac3)
phibac4<-merge(phibac3,phibac2,"otu2")

phibac4[order(phibac4$pi,decreasing=T),]








###### Abundance of hub/connector species ######
#bactabun1b from the Plotting_Bargraph_Regression.R file has all species relative abundances (bact,its,18S). made up of the below dataframes:
datBacr3fotu3[1:10,30:40]
datEukN5fotu3
datEukS4fotu3
datITS3fotu3[1:10,30:40]

hubslo
hubsme
hubshi

zipioutputlo
zipioutputme
zipioutputhi

#mean relative abundances of all bacteria in dataset entering network analysis
relabunlomehib<-aggregate.data.frame(datBacr3fotu3[,32:3803],by=list(datBacr3fotu3$lomehi),mean)
rownames(relabunlomehib)<-relabunlomehib$Group.1
relabunlomehib<-relabunlomehib[,-1]
relabunlomehib<-relabunlomehib*100
relabunlomehib[,1:40]
relabunlomehibt<-t(relabunlomehib)
relabunlomehibt[1:40,]
relabunlomehibt2<-data.frame(otu=rownames(relabunlomehibt),relabunlomehibt)
head(relabunlomehibt2)

zipioutputlorelabun<-merge(zipioutputlo,relabunlomehibt2) #this only selects the bacteria (the intersection of otu in both data sets)
head(zipioutputlorelabun)
zipioutputmerelabun<-merge(zipioutputme,relabunlomehibt2) #this only selects the bacteria (the intersection of otu in both data sets)
head(zipioutputlorelabun)
zipioutputhirelabun<-merge(zipioutputhi,relabunlomehibt2) #this only selects the bacteria (the intersection of otu in both data sets)
head(zipioutputlorelabun)

plot(zipioutputlorelabun$lo,zipioutputlorelabun$zi,pch=1,bg=1)
points(zipioutputmerelabun$me,zipioutputmerelabun$zi,pch=1,bg=1)
points(zipioutputhirelabun$hi,zipioutputhirelabun$zi,pch=1,bg=1)
abline(h=2.5)

#organizing dataframe for ggplot
zipioutputlorelabun2<-zipioutputlorelabun[,-c(6,8)];colnames(zipioutputlorelabun2)[6]<-"relabun"
zipioutputmerelabun2<-zipioutputmerelabun[,-c(6,7)];colnames(zipioutputmerelabun2)[6]<-"relabun"
zipioutputhirelabun2<-zipioutputhirelabun[,-c(7,8)];colnames(zipioutputhirelabun2)[6]<-"relabun"

zipioutputallrelabun2<-rbind(zipioutputlorelabun2,zipioutputmerelabun2,zipioutputhirelabun2)
zipioutputallrelabun2$density<-factor(zipioutputallrelabun2$density,levels=c("Low",'Medium','High'))
zipioutputallrelabun3 <- zipioutputallrelabun2 %>%
  gather(zipi,value,zi:pi)
zipioutputallrelabun3$zipi<-factor(zipioutputallrelabun3$zipi,levels=c("zi","pi"))
zipilines<-data.frame(zipi = c("zi", "pi"), Z = c(2.5, .62))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/zipirelabun.pdf",width=6.4, height=3) #for two panels width=7, height=3.5)
ggplot(zipioutputallrelabun3,aes(x=relabun,y=value,color=density))+# as.numeric(fert),color=species
  labs(x="Relative abundance (%)",y="Within-module connectivity (Zi)")+
  theme_classic()+
  #scale_colour_grey()+
  scale_color_manual(values=c(gray(.7), gray(.5),gray(.3)))+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  #xlim(0,1)+
  #ylim(-2,6)+
  geom_point(size=1.4)+
  geom_hline(data=zipilines,aes(yintercept=Z))+#geom_hline(yintercept=c(2.5,.62))
  facet_wrap(~zipi,scales="free")
dev.off()

zipioutputallrelabun2hubs<-subset(zipioutputallrelabun2,zipioutputallrelabun2$zi>2.5)
min(zipioutputallrelabun2hubs$relabun)
max(zipioutputallrelabun2hubs$relabun)

zipioutputallrelabun2cons<-subset(zipioutputallrelabun2,zipioutputallrelabun2$pi>.62)
min(zipioutputallrelabun2cons$relabun)
max(zipioutputallrelabun2cons$relabun)




#Fungi in hubs/connectors
datITS3fotu3
#mean relative abundances of all fungi in dataset entering network analysis
relabunlomehibf<-aggregate.data.frame(datITS3fotu3[,32:1156],by=list(datITS3fotu3$lomehi),mean)
rownames(relabunlomehibf)<-relabunlomehibf$Group.1
relabunlomehibf<-relabunlomehibf[,-1]
relabunlomehibf<-relabunlomehibf*100
relabunlomehibf[,1:40]
relabunlomehibft<-t(relabunlomehibf)
relabunlomehibft[1:40,]
relabunlomehibft2<-data.frame(otu=rownames(relabunlomehibft),relabunlomehibft)
head(relabunlomehibft2)

zipioutputlorelabunf<-merge(zipioutputlo,relabunlomehibft2) #this only selects the fungi (the intersection of otu in both data sets)
head(zipioutputlorelabunf)
zipioutputmerelabunf<-merge(zipioutputme,relabunlomehibft2) #this only selects the fungi (the intersection of otu in both data sets)
head(zipioutputmerelabunf)
zipioutputhirelabunf<-merge(zipioutputhi,relabunlomehibft2) #this only selects the fungi (the intersection of otu in both data sets)
head(zipioutputhirelabunf)

plot(zipioutputlorelabunf$lo,zipioutputlorelabunf$zi,pch=1,bg=1)
points(zipioutputmerelabunf$me,zipioutputmerelabunf$zi,pch=1,bg=1)
points(zipioutputhirelabunf$hi,zipioutputhirelabunf$zi,pch=1,bg=1)
abline(h=2.5)

plot(zipioutputlorelabunf$lo,zipioutputlorelabunf$pi,pch=1,bg=1,xlim=c(0,13),ylim=c(0,1))
points(zipioutputmerelabunf$me,zipioutputmerelabunf$pi,pch=1,bg=1)
points(zipioutputhirelabunf$hi,zipioutputhirelabunf$pi,pch=1,bg=1)
abline(h=.62)

#organizing dataframe for ggplot
zipioutputlorelabunf2<-zipioutputlorelabunf[,-c(6,8)];colnames(zipioutputlorelabunf2)[6]<-"relabun"
zipioutputmerelabunf2<-zipioutputmerelabunf[,-c(6,7)];colnames(zipioutputmerelabunf2)[6]<-"relabun"
zipioutputhirelabunf2<-zipioutputhirelabunf[,-c(7,8)];colnames(zipioutputhirelabunf2)[6]<-"relabun"

zipioutputallrelabunf2<-rbind(zipioutputlorelabunf2,zipioutputmerelabunf2,zipioutputhirelabunf2)
zipioutputallrelabunf2$density<-factor(zipioutputallrelabunf2$density,levels=c("Low",'Medium','High'))
zipioutputallrelabunf3 <- zipioutputallrelabunf2 %>%
  gather(zipi,value,zi:pi)
zipioutputallrelabunf3$zipi<-factor(zipioutputallrelabunf3$zipi,levels=c("zi","pi"))

min(zipioutputallrelabunf3$relabun)
max(zipioutputallrelabunf3$relabun)









###### Number of plants bacteria and euks in networks and photosynthetic/nonphotosynthetic members  #####
verticesgraphlo<-as.data.frame(rownames(as.matrix(V(graphlo))))
verticesgraphlo2<-labelsall[which(labelsall$otu%in%verticesgraphlo[,1]),]
length(which(verticesgraphlo2$group=="NonphotosyntheticBacteria"|verticesgraphlo2$group=="PhotosyntheticBacteria"))
length(which(verticesgraphlo2$group=="Fungi"))
length(which(verticesgraphlo2$group=="PhotosyntheticEukaryota"|verticesgraphlo2$group=="NonphotosyntheticEukaryota"))
length(which(verticesgraphlo2$group=="Plant"))
length(which(verticesgraphlo2$group=="PhotosyntheticEukaryota"|verticesgraphlo2$group=="PhotosyntheticBacteria")) #27.8% are photosynthetic
length(which(verticesgraphlo2$group=="NonphotosyntheticBacteria"|verticesgraphlo2$group=="NonphotosyntheticEukaryota"|verticesgraphlo2$group=="Fungi")) 

verticesgraphme<-as.data.frame(rownames(as.matrix(V(graphme))))
verticesgraphme2<-labelsall[which(labelsall$otu%in%verticesgraphme[,1]),]
length(which(verticesgraphme2$group=="NonphotosyntheticBacteria"|verticesgraphme2$group=="PhotosyntheticBacteria"))
length(which(verticesgraphme2$group=="Fungi"))
length(which(verticesgraphme2$group=="PhotosyntheticEukaryota"|verticesgraphme2$group=="NonphotosyntheticEukaryota"))
length(which(verticesgraphme2$group=="Plant"))
length(which(verticesgraphme2$group=="PhotosyntheticEukaryota"|verticesgraphme2$group=="PhotosyntheticBacteria"|verticesgraphme2$group=="Plant")) #30.6% are photosynthetic
length(which(verticesgraphme2$group=="NonphotosyntheticBacteria"|verticesgraphme2$group=="NonphotosyntheticEukaryota"|verticesgraphme2$group=="Fungi"))

verticesgraphhi<-as.data.frame(rownames(as.matrix(V(graphhi))))
verticesgraphhi2<-labelsall[which(labelsall$otu%in%verticesgraphhi[,1]),]
length(which(verticesgraphhi2$group=="NonphotosyntheticBacteria"|verticesgraphhi2$group=="PhotosyntheticBacteria"))
length(which(verticesgraphhi2$group=="Fungi"))
length(which(verticesgraphhi2$group=="PhotosyntheticEukaryota"|verticesgraphhi2$group=="NonphotosyntheticEukaryota"))
length(which(verticesgraphhi2$group=="Plant"))
length(which(verticesgraphhi2$group=="Metazoa"))
length(which(verticesgraphhi2$group=="PhotosyntheticEukaryota"|verticesgraphhi2$group=="PhotosyntheticBacteria"|verticesgraphhi2$group=="Plant")) #21.4 are photosynthetic
length(which(verticesgraphhi2$group=="NonphotosyntheticBacteria"|verticesgraphhi2$group=="NonphotosyntheticEukaryota"|verticesgraphhi2$group=="Fungi"|verticesgraphhi2$group=="Metazoa"))


#to get proportion of photosynthetic taxa in the input dataset, I would need to get a list of otu names and then merge that list with the labelfile- bactabun4 has all the microbes that went into the network (and gen/spec analysis which is the same cutoff)
pnonp<-data.frame(otu=rownames(bactabun4))
pnonp2<-merge(pnonp,labelsall,all.x=F,all.y=F,sort=F)
length(which(pnonp2$group%in%c("PhotosyntheticBacteria","PhotosyntheticEukaryota","Plant")))
length(which(pnonp2$group%in%c("Fungi","Metazoa","NonphotosyntheticBacteria","NonphotosyntheticEukaryota")))#,"unknown"
1008/(1008+5611)
#Then I guess I would have to take a random sample of 600 (or mean of taxa in 3 networks), calculate percent photosynthetic, do this 1000 times to get 95% CI, then see if those differences are significant

#this is on the whole data set, not the input datafile
length(which(labelsall$group%in%c("PhotosyntheticBacteria","PhotosyntheticEukaryota","Plant")))
length(which(labelsall$group%in%c("Fungi","Metazoa","NonphotosyntheticBacteria","NonphotosyntheticEukaryota")))#,"unknown"
4520/(4520+30586)



#number of bacteria, fungi, euks in networks as percentage of input
#bacteria
mean(c(546,551,651))/3772

#fungi
mean(c(35,29,25))/1125

#small euks
mean(c(40,33,34))/1498

#large euks - mesofauna
mean(c(0,0,3))/233


#number of bacteria fungi euks in networks for lo me hi separately, as percentage of input of lo me hi separately
datBacr3fotu3
datEukN5fotu3
datEukS4fotu3
datITS3fotu3
plantcomp2

datBacr3fotu3r<-datBacr3fotu3
datBacr3fotu3r$lomehi<-factor(datBacr3fotu3r$lomehi,levels=c("lo","me","hi"))
datBacr3fotu3r2<-aggregate.data.frame(datBacr3fotu3r[,32:3803], by=list(datBacr3fotu3r$lomehi),sum)
datBacr3fotu3r2[,1:10]
temp<-data.frame(lomehi=c("lo","me","hi"),richness=rowSums(datBacr3fotu3r2[,2:3773]>0))
c(546,551,651)/temp$richness

# datBacr3fotu3r2<-ifelse(datBacr3fotu3r[,2:3773]>0,1,0) #check that it is correct
# rowSums(datBacr3fotu3r2)

datITS3fotu3r<-datITS3fotu3
datITS3fotu3r$lomehi<-factor(datITS3fotu3r$lomehi,levels=c("lo","me","hi"))
datITS3fotu3r2<-aggregate.data.frame(datITS3fotu3r[,32:1156], by=list(datITS3fotu3r$lomehi),sum)
datITS3fotu3r2[,1:10]
temp<-data.frame(lomehi=c("lo","me","hi"),richness=rowSums(datITS3fotu3r2[,2:1126]>0))
c(35,29,25)/temp$richness

datEukS4fotu3
datEukS4fotu3r<-datEukS4fotu3
datEukS4fotu3r$lomehi<-factor(datEukS4fotu3r$lomehi,levels=c("lo","me","hi"))
datEukS4fotu3r2<-aggregate.data.frame(datEukS4fotu3r[,32:1520], by=list(datEukS4fotu3r$lomehi),sum)
datEukS4fotu3r2[,1:10]
temp<-data.frame(lomehi=c("lo","me","hi"),richness=rowSums(datEukS4fotu3r2[,2:1490]>0))
c(40,33,34)/temp$richness

datEukN5fotu3
datEukN5fotu3r<-datEukN5fotu3
datEukN5fotu3r$lomehi<-factor(datEukN5fotu3r$lomehi,levels=c("lo","me","hi"))
datEukN5fotu3r2<-aggregate.data.frame(datEukN5fotu3r[,32:264], by=list(datEukN5fotu3r$lomehi),sum)
datEukN5fotu3r2[,1:10]
temp<-data.frame(lomehi=c("lo","me","hi"),richness=rowSums(datEukN5fotu3r2[,2:234]>0))
c(0,0,3)/temp$richness

plantcomp2
plantcomp2r<-comm.dataALL[,6651:6704]
plantcomp2r2<-plantcomp2r[,colSums(plantcomp2r>0)>2]
plantcomp2r3<-data.frame(comm.dataALL[,1:31],plantcomp2r2)
plantcomp2r3$lomehi<-factor(plantcomp2r3$lomehi,levels=c("lo","me","hi"))
plantcomp2r4<-aggregate.data.frame(plantcomp2r3[,32:82], by=list(plantcomp2r3$lomehi),sum)
plantcomp2r4[,1:10]
temp<-data.frame(lomehi=c("lo","me","hi"),richness=rowSums(plantcomp2r4[,2:52]>0))
c(0,2,3)/temp$richness



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
