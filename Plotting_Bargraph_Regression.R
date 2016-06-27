
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
  mutate(type="D. Large Eukaryotes")

names(which(colSums(dat99Met2fNem2k2[,32:39])>.7))
relNem<-dat99Met2fNem2k2 %>% 
  dplyr::select(lomehi,Sample_name, Plant_Div, Plant_Dens,Bacterial_feeder,Fungal_feeder,Omnivore,Plant_parasite,Root_associate) %>%
  gather(Taxa,abun,Bacterial_feeder,Fungal_feeder,Omnivore,Plant_parasite,Root_associate) %>%
  filter(is.na(abun)==F) %>%
  mutate(type="E. Nematodes")

relALL<-rbind(relBac,relITS,relEukS,relEukN,relNem)#
head(relALL)

plotdata<-relALL %>%
  mutate(typeTaxa=paste(type,Taxa)) %>%
  group_by(Taxa,lomehi,type,typeTaxa) %>%
  summarise(mean_abun = mean(abun),se_abun=std.error(abun)) 
  #%>%filter(mean_abun>.04)
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
  facet_wrap(~type,nrow=3,scales="free")
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


#doing anova on glomeromycota
options(contrasts=c("contr.helmert","contr.poly"));options("contrasts")

glom<-relALL %>%
  filter(Taxa=="Glomeromycota")
anova(lm(abun~lomehi,data=glom))






