
biogeo<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Biogeochemistry/CN_enzymes_pH_moisture_whc_Niwot2015t.csv")
head(biogeo)
#note on dataset - I updated this because for the data from microbial biomass (gravimetric moisture, IN, DOC, etc), sample 33 and 34 were mixed up. origianlly sample 33 was missing from the first dataset dorota gave me (CN_enzymes_pH_moisture_whc_Niwot2015.xlsx), and the nubmers fom 34 were really the values for 33. Then I got the All Data file from Dorota and this cleared it up b/c it had both samples 33 adn 34 in it, I also double checked with the Biogeochemistry_Niwot_2015.xlsx file where the calculations were done.

#dataframe that has lomehi in it, datBacr3fotu

biogeo2<-merge(biogeo,datBacr3fotu[,1:31],by="X.SampleID")

dim(biogeo2)

#I could replace all negative numbers with 0, not sure if I should do this, since they are relative, but it doesn't really make sense to have negative microbial biomass or inorganic N
#biogeo2[biogeo2<0]<-0

biogeom<-biogeo2%>%
  select(pH.x:TC,lomehi)%>%
  group_by(lomehi)%>%
  #summarise_each(funs(mean(., na.rm=T),std.error),pH.x:MicN)%>%
  #summarise_each(funs(max(.,na.rm=T)-min(.,na.rm=T)),pH.x:MicN)%>%
  summarise_at(vars(pH.x:TC),funs(sd(.,na.rm=T)/mean(.,na.rm=T)))%>%
  mutate(lomehi=factor(lomehi,levels=c("lo","me","hi")))%>%
  gather(variable,range,pH.x:TC)
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



biogeom<-biogeo2%>%
  select(pH.x:TC,lomehi)%>%
  group_by(lomehi)%>%
  #summarise_each(funs(mean(., na.rm=T),std.error),pH.x:MicN)%>%
  #summarise_each(funs(max(.,na.rm=T)-min(.,na.rm=T)),pH.x:MicN)%>%
  summarise_at(vars(pH.x:TC),funs(mean(.,na.rm=T)))%>%
  mutate(lomehi=factor(lomehi,levels=c("lo","me","hi")))%>%
  gather(variable,range,pH.x:TC)
as.data.frame(biogeom)

plot(biogeo2$Plant_Dens,biogeo2$TC)
summary(lm(TC~Plant_Dens,data=biogeo2))
library(Hmisc)
rcorr(cbind(Plant_Den=biogeo2$Plant_Den,TC=biogeo2$TC,TN=biogeo2$TN))
plot(biogeo2$Plant_Div,biogeo2$TC)

plot(biogeo2$Plant_Dens,biogeo2$TN)
summary(lm(TN~Plant_Dens,data=biogeo2))
plot(biogeo2$Plant_Div,biogeo2$TN)
