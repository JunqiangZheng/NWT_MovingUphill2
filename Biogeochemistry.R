
biogeo<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Biogeochemistry/CN_enzymes_pH_moisture_whc_Niwot2015t.csv")
head(biogeo)

#dataframe that has lomehi in it, datBacr3fotu

biogeo2<-merge(biogeo,datBacr3fotu[,1:31],by="X.SampleID")

dim(biogeo2)

biogeo2[biogeo2<0]<-0

biogeom<-biogeo2%>%
  select(pH.x:MicN,lomehi)%>%
  group_by(lomehi)%>%
  #summarise_each(funs(mean(., na.rm=T),std.error),pH.x:MicN)%>%
  #summarise_each(funs(max(.,na.rm=T)-min(.,na.rm=T)),pH.x:MicN)%>%
  summarise_each(funs(sd(.,na.rm=T)/mean(.,na.rm=T)),pH.x:MicN)%>%
  mutate(lomehi=factor(lomehi,levels=c("lo","me","hi")))%>%
  gather(variable,range,pH.x:MicN)
as.data.frame(biogeom)

range(biogeo2$pH.x,na.rm=T)

ggplot(biogeom,aes(x=lomehi,y=range,color=variable))+
  labs(x = "",y="Range")+
  geom_line(stat = "identity", position = "identity",size=.5)+
  geom_point(size=2)+
  facet_wrap(~variable,scales="free")

none of the enzymes have highest variability in high plants
Highest in high plants: DOC, IN (barely), MicN, NO3, TDN, WHC
Highest somewhere else: MicC, moisture, NH4, pH

