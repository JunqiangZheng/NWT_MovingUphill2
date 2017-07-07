
###### Alpha diversity ######

#pd16S<-pd(as.matrix(datBacr3fotu[,-c(1:31)]),phy_tree(datBac3f),include.root=TRUE)
pdBac$X.SampleID<-rownames(pdBac)
head(datBacr3fotu)[,1:34]
pdBac2<-merge(pdBac,datBacr3fotu[,1:31],"X.SampleID")
pdBac2$type<-"1Bacteria"

#richITS<-as.data.frame(rowSums(datITS3fotu[,-c(1:31)]>0))
colnames(richITS)<-"PD"
richITS$SR<-richITS$PD
richITS$X.SampleID<-rownames(richITS)
richITS2<-merge(richITS,datITS3fotu[,1:31],"X.SampleID")
richITS2$type<-"2Fungi"

#pdEukN<-pd(as.matrix(datEukN5fotu[,-c(1:31)]),phy_tree(datEukN5f),include.root=TRUE) 
pdEukN$X.SampleID<-rownames(pdEukN)
head(datEukN5fotu)[,1:34]
pdEukN2<-merge(pdEukN,datEukN5fotu[,1:31],"X.SampleID")
pdEukN2$type<-"4Large Eukaryotes"

#pdEukS<-pd(as.matrix(datEukS4fotu[,-c(1:31)]),phy_tree(datEukS4f),include.root=TRUE) 
pdEukS$X.SampleID<-rownames(pdEukS)
head(datEukS4fotu)[,1:34]
pdEukS2<-merge(pdEukS,datEukS4fotu[,1:31],"X.SampleID")
pdEukS2$type<-"3Small Eukaryotes"

#pdMet<-pd(as.matrix(dat99Met2fotu[,-c(1:31)]),phy_tree(dat99Met2f),include.root=TRUE) 
pdMet$X.SampleID<-rownames(pdMet)
head(dat99Met2fotu)[,1:34]
pdMet2<-merge(pdMet,dat99Met2fotu[,1:31],"X.SampleID")
pdMet2$type<-"Metazoa"

#pdNem<-pd(as.matrix(dat99Met2fNemotu[,-c(1:31)]),phy_tree(dat99Met2fNem),include.root=TRUE)
pdNem$X.SampleID<-rownames(pdNem)
head(dat99Met2fNemotu)[,1:34]
pdNem2<-merge(pdNem,dat99Met2fNemotu[,1:31],"X.SampleID")
pdNem2$type<-"5Nematodes"

plot(pdNem2$Plant_Dens,pdNem2$PD)

richdata<-rbind(pdBac2,richITS2,pdEukN2,pdEukS2,pdNem2)

#fig with euks and bac
pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/bacfunsleuknemPDvsplantdensity.pdf",width=7, height=4.7) #for two panels width=7, height=3.5)
ggplot(richdata,aes(x=log10(Plant_Dens+1),y=PD))+# as.numeric(fert),color=species
  labs(x="Plant density",y="Diversity")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=12),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.5))+
  geom_point(size=1.4)+
  geom_smooth(method=lm,se=F,size=.8,color="black") +
  #geom_smooth(method=lm,se=F,size=.8,color="black",formula = y ~ poly(x, 2)) +
  facet_wrap(~type,scales="free")
dev.off()

summary(lm(pdBac2$PD~log((pdBac2$Plant_Dens+1),base=10)))
summary(lm(richITS2$PD~log((richITS2$Plant_Dens+1),base=10)))
summary(lm(pdEukS2$PD~log((pdEukS2$Plant_Dens+1),base=10)))
summary(lm(pdEukN2$PD~log((pdEukN2$Plant_Dens+1),base=10)))
summary(lm(pdNem2$PD~log((pdNem2$Plant_Dens+1),base=10)))

#Fitting curves
b1<-lm(pdBac2$PD~log((pdBac2$Plant_Dens+1),base=10))
b2<-nls(PD~a+b*log((Plant_Dens+1),base=10)+c*log((Plant_Dens+1),base=10)^2,start=c(a=1,b=1,c=1),data=pdBac2)

n1<-lm(pdNem2$PD~log((pdNem2$Plant_Dens+1),base=10))
n2<-nls(PD~a+b*log((Plant_Dens+1),base=10)+c*log((Plant_Dens+1),base=10)^2,start=c(a=1,b=1,c=1),data=pdNem2)
AIC(n1)
AIC(n2)








###### Change in relative abundance ######

colnames(datBac3fk2)[78]<-"WCHB1.60"
names(which(colSums(datBac3fk2[,32:78])>4))
relBac<-datBac3fk2 %>% 
  dplyr::select(lomehi,Sample_name, Plant_Div, Plant_Dens,Acidobacteria,Actinobacteria,Bacteroidetes,Chloroflexi,Cyanobacteria,Planctomycetes,Proteobacteria,Verrucomicrobia) %>%
  gather(Taxa,abun,Acidobacteria:Verrucomicrobia) %>%
  mutate(type="A. Bacteria")

names(which(colSums(datITS3fk2[,32:39])>1))
relITS<-datITS3fk2 %>% 
  dplyr::select(lomehi,Sample_name, Plant_Div, Plant_Dens,Ascomycota,Basidiomycota,Glomeromycota,Zygomycota) %>%
  gather(Taxa,abun,Ascomycota,Basidiomycota,Glomeromycota,Zygomycota) %>%
  mutate(type="B. Fungi")

names(which(colSums(datEukS4fk2[,32:46])>4))#,Nonphotosynthetic_Excavata,Photosynthetic_Stramenopiles,
relEukS<-datEukS4fk2 %>% 
  dplyr::select(lomehi,Sample_name, Plant_Div, Plant_Dens,Amoebozoa,Archaeplastida,Nonphotosynthetic_Alveolata,Rhizaria) %>%
  gather(Taxa,abun,Amoebozoa,Archaeplastida,Nonphotosynthetic_Alveolata,Rhizaria) %>%
  mutate(type="C. Small Eukaryotes")

names(which(colSums(datEukN5fk2[,32:43])>.7))
relEukN<-datEukN5fk2 %>% 
  dplyr::select(lomehi,Sample_name, Plant_Div, Plant_Dens,Arthropoda,Nematoda,Rotifera,Tardigrada) %>%
  gather(Taxa,abun,Arthropoda,Nematoda,Rotifera,Tardigrada) %>%
  mutate(type="D. Soil mesofauna")

names(which(colSums(dat99Met2fNem2k2[,32:39])>.7))
relNem<-dat99Met2fNem2k2 %>% 
  dplyr::select(lomehi,Sample_name, Plant_Div, Plant_Dens,Bacterial_feeder,Fungal_feeder,Omnivore,Plant_parasite,Root_associate) %>%
  gather(Taxa,abun,Bacterial_feeder,Fungal_feeder,Omnivore,Plant_parasite,Root_associate) %>%
  filter(is.na(abun)==F) %>%
  mutate(type="E. Nematodes")

relALL<-rbind(relBac,relITS,relEukS,relEukN,relNem)#
head(relALL)

#plotdata<-relALL %>%
#  mutate(typeTaxa=paste(type,Taxa)) %>%
#  group_by(Taxa,lomehi,type,typeTaxa) %>%
#  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
#  #%>%filter(mean_abun>.04)

#this was weird, maybe something changed in ggplot or dplyr because the colors were messing up and it was listing the legend in alfabetical order by taxa rather than the order in the "plotdata" dataframe. the workaroudn was to set the levels of plotdata$Taxa so they were correct
plotdata<-relALL %>%
  mutate(typeTaxa=paste(type,Taxa)) %>%
  group_by(typeTaxa,Taxa,lomehi,type) %>%
  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
plotdata$Taxa<-factor(plotdata$Taxa,levels=unique(plotdata$Taxa))

as.data.frame(plotdata)
plotdata$lomehi<-factor(plotdata$lomehi,levels=c("lo","me","hi"))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/relabuntaxavsplantdensitygroupsBFSLEN.pdf",width=6.5,height=6)#,width=4.3, height=5.3
ggplot(plotdata,aes(x=lomehi,y=mean_abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.15,size=.5)+
  scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()

#8 bacteria, 4 fungi, 4 small euks, 4 large euks, 5 nem
mycols<-c("#4BC366",
          "#D9A125",
          "#659125",
          "#6768A3",
          "#5C426C",
          "#D185E0",
          "#6F94DE",
          "#B4405E",
          "#D9A125",#yellow
          "#B4405E",#red
          "#659125",
          "#6768A3",
          "#5C426C",
          "#D185E0",
          "#6F94DE",
          "#4BC366",
          "#4BC366", #light green
          "#B4405E", #red
          "#5C426C", #dark purple
          "#D9A125", #yellow
          "#659125", #darker green
          "#D185E0", #light purple
          "#6768A3", #medium purple
          "#D9A125", #yellow
          "#6F94DE") #blue)

#scatter plots
head(relALL)

plotdata<-relALL %>%
  mutate(typeTaxa=paste(type,Taxa))
plotdata$Taxa<-factor(plotdata$Taxa,levels=unique(plotdata$Taxa))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/relabuntaxavsplantdensitygroupsBFSLENscatter.pdf",width=6.5,height=6)#,width=4.3, height=5.3
ggplot(plotdata,aes(x=log10(Plant_Dens+1),y=abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  scale_color_manual(values=mycols) +
  geom_point(size=.2)+
  geom_smooth(method=lm,se=F,size=.8) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()


#Doing anova on all of the above taxa groups
length(unique(relALL$Taxa))
anovaoutput<-data.frame(Taxa=rep(NA,25),F=rep(NA,25),P=rep(NA,25))
for(i in 1:25){
  current.taxa<-unique(relALL$Taxa)[i]
  temp<-relALL %>%
    filter(Taxa==current.taxa)
  mod<-anova(lm(abun~lomehi,data=temp))
  anovaoutput[i,1]<-current.taxa
  anovaoutput[i,2]<-mod$`F value`[1]
  anovaoutput[i,3]<-mod$`Pr(>F)`[1]
}
anovaoutput$qval<-p.adjust(anovaoutput$P,method="fdr")
anovaoutput$qval<-format(anovaoutput$qval,scientific=F)

#doing anova on glomeromycota
options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")

glom<-relALL %>%
  filter(Taxa=="Glomeromycota")
anova(lm(abun~lomehi,data=glom))

chyt<-datITS3fk2 %>% 
  dplyr::select(lomehi,Sample_name, Plant_Div, Plant_Dens,Chytridiomycota) %>%
  gather(Taxa,abun,Chytridiomycota) %>%
  mutate(type="B. Fungi")%>%
  group_by(lomehi)%>%
  summarise(mean_abun=mean(abun),se_abun=std.error(abun))
chyt$lomehi<-factor(chyt$lomehi, levels=c("lo","me","hi"))

ggplot(chyt,aes(x=lomehi,y=mean_abun))+
  labs(x = "",y="Relative abundance")+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.15,size=.5)
    
    
#looking into Ktedonobacteria in Chloroflexi
#I would have to regenerate the dataset at the class level (above datasets are at the phylum level)







##### Specialists/generalists ######
#input data, note these are not in order by sample
datBacr3fotu
datEukN5fotu
datEukS4fotu
datITS3fotu

#Bacteria specialists/generalists
bactspgenlo1<-subset(datBacr3fotu,lomehi=="lo")
bactspgenlo2<-bactspgenlo1[,-c(1:31)]
bactspgenlo3<-bactspgenlo2[,which(colSums(bactspgenlo2)>0)]
bactspgenlo<-colnames(bactspgenlo3)
length(bactspgenlo) #12023

bactspgenme1<-subset(datBacr3fotu,lomehi=="me")
bactspgenme2<-bactspgenme1[,-c(1:31)]
bactspgenme3<-bactspgenme2[,which(colSums(bactspgenme2)>0)]
bactspgenme<-colnames(bactspgenme3)
length(bactspgenme) #14024

bactspgenhi1<-subset(datBacr3fotu,lomehi=="hi")
bactspgenhi2<-bactspgenhi1[,-c(1:31)]
bactspgenhi3<-bactspgenhi2[,which(colSums(bactspgenhi2)>0)]
bactspgenhi<-colnames(bactspgenhi3)
length(bactspgenhi)#16688


#generalist bacteria

temp<-intersect(bactspgenlo,bactspgenme)
bacttaxagen<-intersect(temp,bactspgenhi)
length(bacttaxagen)
#5911 taxa are shared across all communities
3*5911/(12023+14024+16699)

temp<-datBacr3fotu %>%
  select(one_of(c("lomehi",bacttaxagen)))%>%
  mutate(sum=rowSums(.[2:5912]))%>%
  group_by(lomehi)%>%
  summarise(mean=mean(sum),se=std.error(sum))%>%
  mutate(lomehi=factor(lomehi,levels=c("lo","me","hi")))

ggplot(temp,aes(x=lomehi,y=mean))+
  labs(x = "",y="Relative abundance")+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.5)


#specialist bacteria

temp<-setdiff(bactspgenlo,bactspgenme)
bacttaxasplo<-setdiff(temp,bactspgenhi)
length(bacttaxasplo)

bactspgenlo3%>%
  select(one_of(bacttaxasplo))%>%
  mutate(sum=rowSums(.))%>%
  select(sum)%>%
  summarise(mean=mean(sum),se=std.error(sum))
  
temp<-setdiff(bactspgenhi,bactspgenme)
bacttaxasphi<-setdiff(temp,bactspgenlo)
length(bacttaxasphi)

bactspgenhi3%>%
  select(one_of(bacttaxasphi))%>%
  mutate(sum=rowSums(.))%>%
  select(sum)%>%
  summarise(mean=mean(sum),se=std.error(sum))

temp<-bactspgenlo3
temp[which(temp==0,arr.ind=T)]<-NA
min(temp,na.rm=T)

library(scatterplot3d) 


#investigating how abundance and frequency plays out across the data set - most species are rare, common species are found in all lo me and hi
bactabun1<-cbind(lomehi=datBacr3fotu$lomehi,datBacr3fotu[,-c(1:31)])
bactabun1[1:10,1:10]
bactabun2<-aggregate.data.frame(bactabun1[,2:25412],by=list(bactabun1$lomehi),sum)
dim(bactabun2)
bactabun2[,1:10]
bactabun3<-t(bactabun2[,-c(1)])
bactabun3[1:10,];colnames(bactabun3)<-c("hi","lo","me")
bactabun3[1:10,]
bactabun4<-as.data.frame(bactabun3)
#scatterplot3d(log(bactabun4$lo+.1),log(bactabun4$me+.1),log(bactabun4$hi+.1),type="h")
head(bactabun4)
bactabun4[ind,]*100

bactfreq1<-data.frame(lomehi=as.character(datBacr3fotu$lomehi),ifelse(datBacr3fotu[,-c(1:31)]>0,1,0))
bactfreq1[1:10,1:10]
bactfreq2<-aggregate.data.frame(bactfreq1[,2:25412],by=list(bactfreq1$lomehi),sum)
dim(bactfreq2)
bactfreq2[,1:10]
bactfreq3<-t(bactfreq2[,-c(1)])
bactfreq3[1:10,];colnames(bactfreq3)<-c("hi","lo","me")
bactfreq3[1:10,]
bactfreq4<-as.data.frame(bactfreq3)
#scatterplot3d(log(bactfreq4$lo+1),log(bactfreq4$me+1),log(bactfreq4$hi+1),type="h")
head(bactfreq4)
ind<-which(rowSums(bactfreq4)>9)
length(ind)
bactfreq4[ind,]

#making the plot from barbaran et al, abun vs. freq
plot(rowSums(bactfreq4),rowSums(bactabun4),log="y")
plot(rowSums(bactfreq4),rowSums(bactabun4))


#specialist = at least 75% (or other percent) of the abundance or frequency of a taxon is in one plant density bracket
#generalist = in all plots and not as above

#for abundance
bactabun5<-bactabun4
bactabun5$hip<-bactabun4$hi/rowSums(bactabun4)
bactabun5$lop<-bactabun4$lo/rowSums(bactabun4)
bactabun5$mep<-bactabun4$me/rowSums(bactabun4)
head(bactabun5)

#for frequency
bactabun5<-bactfreq4
bactabun5$hip<-bactfreq4$hi/rowSums(bactfreq4)
bactabun5$lop<-bactfreq4$lo/rowSums(bactfreq4)
bactabun5$mep<-bactfreq4$me/rowSums(bactfreq4)
head(bactabun5)

#then continue here for everything
bactabun5$group<-NA
ind<-which(bactabun5$hip>.75)
bactabun5$group[ind]<-'hi'
ind<-which(bactabun5$lop>.75)
bactabun5$group[ind]<-'lo'
ind<-which(bactabun5$mep>.75)
bactabun5$group[ind]<-'me'

bactabun5$class<-NA
ind<-which(rowSums(bactabun5[,1:3]>0)==3)
bactabun5$class[ind]<-"gen"
ind<-which(bactabun5$group%in%c("lo","me","hi"))
bactabun5$class[ind]<-"spe"
#bactabun5[which(rowSums(bactabun5[,1:3]>0)==3),] #just a check that there are specialists that are present in all three groups

#generalists
bactgen<-rownames(bactabun5[which(bactabun5$class=="gen"),])

temp<-datBacr3fotu %>%
  select(one_of(c("lomehi",bactgen)))%>%
  mutate(sum=rowSums(.[2:length(bactgen)+1]))%>%
  group_by(lomehi)%>%
  summarise(mean=mean(sum),se=std.error(sum))%>%
  mutate(lomehi=factor(lomehi,levels=c("lo","me","hi")))

ggplot(temp,aes(x=lomehi,y=mean))+
  labs(x = "",y="Relative abundance")+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.5)

#specialists
bacttaxasplo<-rownames(bactabun5[which(bactabun5$group=="lo"),])
length(bacttaxasplo)

templo<-datBacr3fotu%>%
  select(one_of(bacttaxasplo))%>%
  mutate(sum=rowSums(.))%>%
  select(sum)%>%
  summarise(mean=mean(sum),se=std.error(sum))

bacttaxaspme<-rownames(bactabun5[which(bactabun5$group=="me"),])
length(bacttaxaspme)

tempme<-datBacr3fotu%>%
  select(one_of(bacttaxaspme))%>%
  mutate(sum=rowSums(.))%>%
  select(sum)%>%
  summarise(mean=mean(sum),se=std.error(sum))

bacttaxasphi<-rownames(bactabun5[which(bactabun5$group=="hi"),])
length(bacttaxasphi)

temphi<-datBacr3fotu%>%
  select(one_of(bacttaxasphi))%>%
  mutate(sum=rowSums(.))%>%
  select(sum)%>%
  summarise(mean=mean(sum),se=std.error(sum))

temp<-rbind(templo,tempme,temphi)
temp$lomehi<-factor(c("lo",'me','hi'),levels=c("lo",'me','hi'))

ggplot(temp,aes(x=lomehi,y=mean))+
  labs(x = "",y="Relative abundance")+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.5)






#what I could do
#take out singletons/rare species
#mess with the 75% limit - down to 60

#take out all doubletons and singletons, and taxa with a summed rel abun <.2% (same filter as for what is going into the networks) - this means use: 
datBacr3fotu3
datEukN5fotu3
datEukS4fotu3
datITS3fotu3

#investigating how abundance and frequency plays out across the data set - most species are rare, common species are found in all lo me and hi
bactabun1<-cbind(lomehi=datBacr3fotu3$lomehi,datBacr3fotu3[,-c(1:31)])
bactabun1<-cbind(lomehi=datEukN5fotu3$lomehi,datEukN5fotu3[,-c(1:31)])
bactabun1<-cbind(lomehi=datEukS4fotu3$lomehi,datEukS4fotu3[,-c(1:31)])
bactabun1<-cbind(lomehi=datITS3fotu3$lomehi,datITS3fotu3[,-c(1:31)])
bactabun1[1:10,1:10]
bactabun2<-aggregate.data.frame(bactabun1[,2:dim(bactabun1)[2]],by=list(bactabun1$lomehi),sum)
dim(bactabun2)
bactabun2[,1:10]
bactabun3<-t(bactabun2[,-c(1)])
bactabun3[1:10,];colnames(bactabun3)<-c("hi","lo","me")
bactabun3[1:10,]
bactabun4<-as.data.frame(bactabun3)
#scatterplot3d(log(bactabun4$lo+.1),log(bactabun4$me+.1),log(bactabun4$hi+.1),type="h")
head(bactabun4)

bactfreq1<-data.frame(lomehi=as.character(datBacr3fotu3$lomehi),ifelse(datBacr3fotu3[,-c(1:31)]>0,1,0))
bactfreq1<-data.frame(lomehi=as.character(datEukN5fotu3$lomehi),ifelse(datEukN5fotu3[,-c(1:31)]>0,1,0))
bactfreq1<-data.frame(lomehi=as.character(datEukS4fotu3$lomehi),ifelse(datEukS4fotu3[,-c(1:31)]>0,1,0))
bactfreq1<-data.frame(lomehi=as.character(datITS3fotu3$lomehi),ifelse(datITS3fotu3[,-c(1:31)]>0,1,0))
bactfreq1[1:10,1:10]
bactfreq2<-aggregate.data.frame(bactfreq1[,2:dim(bactabun1)[2]],by=list(bactfreq1$lomehi),sum)
dim(bactfreq2)
bactfreq2[,1:10]
bactfreq3<-t(bactfreq2[,-c(1)])
bactfreq3[1:10,];colnames(bactfreq3)<-c("hi","lo","me")
bactfreq3[1:10,]
bactfreq4<-as.data.frame(bactfreq3)
#scatterplot3d(log(bactfreq4$lo+1),log(bactfreq4$me+1),log(bactfreq4$hi+1),type="h")
head(bactfreq4)
ind<-which(rowSums(bactfreq4)>9)
length(ind)
bactfreq4[ind,]

#making the plot from barbaran et al, abun vs. freq
plot(rowSums(bactfreq4),rowSums(bactabun4),log="y")
plot(rowSums(bactfreq4),rowSums(bactabun4))



#specialist = at least 60% of the abundance of a taxon is in one plant density bracket
#generalist = in all plots and not as above

bactabun5<-bactabun4
bactabun5$hip<-bactabun4$hi/rowSums(bactabun4)
bactabun5$lop<-bactabun4$lo/rowSums(bactabun4)
bactabun5$mep<-bactabun4$me/rowSums(bactabun4)
head(bactabun5)

#for frequency
bactabun5<-bactfreq4
bactabun5$hip<-bactfreq4$hi/rowSums(bactfreq4)
bactabun5$lop<-bactfreq4$lo/rowSums(bactfreq4)
bactabun5$mep<-bactfreq4$me/rowSums(bactfreq4)
head(bactabun5)

bactabun5$group<-NA
ind<-which(bactabun5$hip>.70)
bactabun5$group[ind]<-'hi'
ind<-which(bactabun5$lop>.70)
bactabun5$group[ind]<-'lo'
ind<-which(bactabun5$mep>.70)
bactabun5$group[ind]<-'me'

bactabun5$class<-NA
ind<-which(rowSums(bactabun5[,1:3]>0)==3)
bactabun5$class[ind]<-"gen"
ind<-which(bactabun5$group%in%c("lo","me","hi"))
bactabun5$class[ind]<-"spe"
#bactabun5[which(rowSums(bactabun5[,1:3]>0)==3),] #just a check that there are specialists that are present in all three groups

bactgen<-rownames(bactabun5[which(bactabun5$class=="gen"),])

tempdata<-datBacr3fotu3
tempdata<-datEukN5fotu3 
tempdata<-datEukS4fotu3 
tempdata<-datITS3fotu3 

temp<-tempdata%>%
  select(one_of(c("lomehi",bactgen)))%>%
  mutate(sum=rowSums(.[2:length(bactgen)+1]))%>%
  group_by(lomehi)%>%
  summarise(mean=mean(sum),se=std.error(sum))%>%
  mutate(lomehi=factor(lomehi,levels=c("lo","me","hi")))

ggplot(temp,aes(x=lomehi,y=mean))+
  labs(x = "",y="Relative abundance")+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.5)

#specialists
bacttaxasplo<-rownames(bactabun5[which(bactabun5$group=="lo"),])
length(bacttaxasplo)

templo<-tempdata%>%
  select(one_of(bacttaxasplo))%>%
  mutate(sum=rowSums(.))%>%
  select(sum)%>%
  summarise(mean=mean(sum),se=std.error(sum))

bacttaxaspme<-rownames(bactabun5[which(bactabun5$group=="me"),])
length(bacttaxaspme)

tempme<-tempdata%>%
  select(one_of(bacttaxaspme))%>%
  mutate(sum=rowSums(.))%>%
  select(sum)%>%
  summarise(mean=mean(sum),se=std.error(sum))

bacttaxasphi<-rownames(bactabun5[which(bactabun5$group=="hi"),])
length(bacttaxasphi)

temphi<-tempdata%>%
  select(one_of(bacttaxasphi))%>%
  mutate(sum=rowSums(.))%>%
  select(sum)%>%
  summarise(mean=mean(sum),se=std.error(sum))

temp<-rbind(templo,tempme,temphi)
temp$lomehi<-factor(c("lo",'me','hi'),levels=c("lo",'me','hi'))

ggplot(temp,aes(x=lomehi,y=mean))+
  labs(x = "",y="Relative abundance")+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean+se, ymin=mean-se),width=.15,size=.5)



#####
plotdata<-relALL %>%
  mutate(typeTaxa=paste(type,Taxa)) %>%
  group_by(typeTaxa,Taxa,lomehi,type) %>%
  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
plotdata$Taxa<-factor(plotdata$Taxa,levels=unique(plotdata$Taxa))

as.data.frame(plotdata)
plotdata$lomehi<-factor(plotdata$lomehi,levels=c("lo","me","hi"))

pdf("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Figs/relabuntaxavsplantdensitygroupsBFSLEN.pdf",width=6.5,height=6)#,width=4.3, height=5.3
ggplot(plotdata,aes(x=lomehi,y=mean_abun,group=typeTaxa,color=Taxa))+
  labs(x = "",y="Relative abundance")+
  theme_classic()+
  theme(line=element_line(size=.3),text=element_text(size=10),strip.background = element_rect(colour="white", fill="white"),axis.line=element_line(color="gray30",size=.3),legend.key.size = unit(.6, "line"))+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  geom_errorbar(aes(ymax = mean_abun+se_abun, ymin=mean_abun-se_abun),width=.15,size=.5)+
  scale_color_manual(values=mycols) +
  facet_wrap(~type,nrow=3,scales="free")+
  guides(col = guide_legend(ncol = 1))
dev.off()


