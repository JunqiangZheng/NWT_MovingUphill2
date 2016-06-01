
#######Read in OTU data#######
#files:
#Euks no metazoa no fungi 97%
#Euks metazoa 97%
#Fungi ITS 97%
#Bacteria 97%

##### Read in euk files, filtered with singletons removed #####
otufileEukS <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_97_S111_OTU_tablefiltsingnonchimericbactarcplantEukSoil2015nometazoanfungising.biom")
otufileEukN <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_97_S111_OTU_tablefiltsingnonchimericbactarcplantEukNematodemetazoancransing.biom")

head(tax_table(otufileEukN))
unique(tax_table(otufileEukS)[,"Rank3"])

#Import mapping and tree file
mapEuk<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/EukBr_Niwot_20072015_All_MapFilenewlomehi.txt")

treeEukS<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_truncate_97_rep_set_filtsingS2015nonchimeras_sinaaln.tre")
treeEukN<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_truncate_97_rep_set_filtsingN2015nonchimeras_sinaaln.tre")

datEukS<-merge_phyloseq(otufileEukS,mapEuk,treeEukS)
datEukN<-merge_phyloseq(otufileEukN,mapEuk,treeEukN)

#remove a few samples, sample 2A for nematodes (rep), sample N.2A has only 3 otus in it and they all are abundant in other samples so it doesn't affect single read count
datEukN2 <- prune_samples(sample_names(datEukN)!="N.2A.2015", datEukN)
#sample 78 for nematodes (it only had 164 reads)
datEukN3 <- prune_samples(sample_names(datEukN2)!="N.78.2015", datEukN2)
#maybe 61 which has 489 reads for the soil dataset
datEukS2 <- prune_samples(sample_names(datEukS)!="S.61.2015", datEukS)

#rarefy and transform to relative abundance
min(sample_sums(datEukS2))#rarefy to 808
min(sample_sums(datEukN3))#rarefy to 692
datEukS3<-rarefy_even_depth(datEukS2,sample.size=min(sample_sums(datEukS2)),rngseed=10,replace=F) #1895 OTUs were removed because they are no longer present in any sample after random subsampling, 
datEukN4<-rarefy_even_depth(datEukN3,sample.size=min(sample_sums(datEukN3)),rngseed=10,replace=F) #2601 OTUs were removed because they are no longer present in any sample after random subsampling, 
datEukS4 = transform_sample_counts(datEukS3, function(x) x/sum(x) )
datEukN5 = transform_sample_counts(datEukN4, function(x) x/sum(x) )
rownames(otu_table(datEukr))[1:10]



#make a label column with kingdom/phylum level labels

#Making groups for labeling graphs
#Fungi (kingdom), 
#Archaeplastida (major_clade), combine the kingdoms Chloroplastida, Rhodophyceae, Haptophyta, 
#Rhizaria (kingdom: unicellular amoeboid euks), 
#Amoebozoa(kingdom),
#Holozoa(kingdom: all animals, not fungi), includes metazoa (animals)
#photosynthetic Alveolata, kingdom (phylum Dinoflagellata: mostly photosynthetic; nonphotosynthetic Protoperidinium, both SL163A10 and Pfiesteria (can if it eats an alga), unknown D244), 
#nonphotosynthetic Alveolata (phyla Ciliophora(predators), protalveolata, apicomplexa (parasites)), BOLA553 is a near apicomplexan, not sure what SCM37C52 is so I put it here
#photosynthetic Excavata major clade (was Discoba kingdom) (class Euglenida: mostly photosynthetic), 
#nonphotosnthetic Excavata (was Discoba) (in discoba, phylum Heterolobosea: parasites, free living, symbiotic, amoeba-like, and phylum Jakoba heterotrophic), superphylum metamonada, doesn't have kingdom, parasites not Ps.
#photosynthetic Cryptophyceae is a kingdom, they have one or two chloroplasts, except for Chilomonas, which has leucoplasts (not pigmented not for Ps) and Goniomonas (formerly Cyathomonas) which lacks plastids entirely.P1-31 (i can't find much info on this, some sites called them picoplankton and most of the cryptophyceae are photosynthetic). Cryptomonadales can't find much about it but assume Ps. __SA1-3C06 is at rank2 but notes in silva say "cryptophyte"
#nonphotosynthetic cryptophyceae - Chilomonas and Goniomonas (formerly Cyathomonas)
#Haptophyta kingdom, I think all are Ps
#Centrohelida - encyclopedia of lif said some species have photosynthetic symbionts but I can't find any info on Ps taxa.Heterophryidae, Mb-5C (can't find any info on, I'll put it with nonPs), M1-18D08(can't find info),Acanthocystidae(can't find info),Pterocystis(can't find info), so I will leave them all in a group called "centrohelida"
#NonPs Stramenopiles: MAST-3, MAST-12, MAST-7, Labyrinthulomycetes,Bicosoecida,Peronosporomycetes,Hyphochytriales,__Incertae_Sedis(were all __Pirsonia which is a parasite of algae)
#Ps Stramenopiles: Pelagophycea, Diatomea, Eustigmatales, Xanthophyceae, Phaeothamniophyceae, Chrysophyceae,Ochromonas,CCI40 (not much info, some places call it is chrysophyte),Synurales, Raphidophyceae, NA1-2A5 (can't find any info so putting it with Ps)

#heterotrophic eukaryota (things without a kingdom): Picozoa (phylum nonphotosynthetic unicellular eukaryotes, only 1 OTU); RT5iin14; Apusomonadidae (nonphotosynthetic protist, unknown kingdom); Kathablepharidae (nonPs protist, unkonwn kingdo);Telonema (nonPs protist genus);Breviatea (amoeba like);DH147-EKD23;__Ancyromonadida;GoC2-B10 (no info);Zeuk77 (no info)

#datEukS4
#datEukN5

head(tax_table(datEukS4))
labelsEukS<-tax_table(datEukS4)[,"Rank3"]#

ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Fungi")
labelsEuk[ind]<-"Fungi"
ind<-which(tax_table(datEukr2)[,"Rank3"]=="__uncultured_Sarcosomataceae")#i'm putting the wierd centrohelid/fungus with fungi
labelsEuk[ind]<-"Fungi"
ind<-which(tax_table(datEukr2)[,"Rank2"]=="__Archaeplastida")
labelsEuk[ind]<-"Archaeplastida"
ind<-which(tax_table(datEukr2)[,"Rank3"]=="__Rhizaria")
labelsEuk[ind]<-"Rhizaria"
ind<-which(tax_table(datEukr2)[,"Rank2"]=="__Amoebozoa")
labelsEuk[ind]<-"Amoebozoa"
ind<-which(tax_table(datEukr2)[,"Rank3"]=="__Holozoa"|tax_table(datEukr2)[,"Rank3"]=="__Metazoa")
labelsEuk[ind]<-"Holozoa"
ind<-which(tax_table(datEukr2)[,"Rank4"]=="__Dinoflagellata")
labelsEuk[ind]<-"Photosynthetic_Alveolata"
ind<-which(tax_table(datEukr2)[,"Rank4"]!="__Dinoflagellata"&tax_table(datEukr2)[,"Rank3"]=="__Alveolata")
labelsEuk[ind]<-"Nonphotosynthetic_Alveolata"
ind<-which(tax_table(datEukr2)[,"Rank3"]=="__Goniomonas")
labelsEuk[ind]<-"Nonphotosynthetic_Chryptophyceae"
ind<-which(tax_table(datEukr2)[,"Rank3"]=="__P1-31"|tax_table(datEukr2)[,"Rank3"]=="__Cryptomonadales"|tax_table(datEukr2)[,"Rank2"]=="__SA1-3C06")
labelsEuk[ind]<-"Photosynthetic_Chryptophyceae"
ind<-which(tax_table(datEukr2)[,"Rank6"]=="__Euglenida")
labelsEuk[ind]<-"Photosynthetic_Excavata"
ind<-which(tax_table(datEukr2)[,"Rank6"]!="__Euglenida"&tax_table(datEukr2)[,"Rank2"]=="__Excavata")
labelsEuk[ind]<-"Nonphotosynthetic_Excavata"
ind<-which(tax_table(datEukr2)[,"Rank2"]=="__Haptophyta")
labelsEuk[ind]<-"Haptophyta"
ind<-which(tax_table(datEukr2)[,"Rank2"]=="__Centrohelida"&tax_table(datEukr2)[,"Rank3"]!="__uncultured_Sarcosomataceae")
labelsEuk[ind]<-"Centrohelida"
ind<-which(tax_table(datEukr2)[,"Rank3"]=="__Stramenopiles"|tax_table(datEukr2)[,"Rank4"]%in%c("__MAST-3","__MAST-12","__MAST-7","__Labyrinthulomycetes","__Bicosoecida","__Peronosporomycetes","__Hyphochytriales,__Incertae_Sedis"))
labelsEuk[ind]<-"Nonphotosynthetic_Stramenopiles"
ind<-which(tax_table(datEukr2)[,"Rank3"]=="__Stramenopiles"|tax_table(datEukr2)[,"Rank4"]%in%c("__Pelagophycea","__Diatomea","__Eustigmatales","__Xanthophyceae","__Phaeothamniophyceae","__Chrysophyceae","__Ochromonas","__CCI40","__Synurales","__Raphidophyceae","__NA1-2A5"))
labelsEuk[ind]<-"Photosynthetic_Stramenopiles"
ind<-which(tax_table(datEukr2)[,"Rank3"]=="__RT5iin14"|tax_table(datEukr2)[,"Rank2"]=="__Picozoa"|tax_table(datEukr2)[,"Rank2"]=="__RT5iin25"|tax_table(datEukr2)[,"Rank3"]=="__Apusomonadidae"|tax_table(datEukr2)[,"Rank2"]=="__Kathablepharidae"|tax_table(datEukr2)[,"Rank3"]=="__Telonema"|tax_table(datEukr2)[,"Rank3"]=="__Breviatea"|tax_table(datEukr2)[,"Rank2"]=="__DH147-EKD23"|tax_table(datEukr2)[,"Rank3"]=="__Ancyromonadida"|tax_table(datEukr2)[,"Rank2"]=="__GoC2-B10"|tax_table(datEukr2)[,"Rank2"]=="__Zeuk77")
labelsEuk[ind]<-"Heterotrophic_Eukarya"
#I could separate excavata into kingdoms and super phyla and also separate archaeplastida into its 3 kingdoms
unique(labelsEuk)
colnames(labelsEuk)<-"labels"

#replace tax table for datEukr2 (relative abun) and datEukr (not relativized)
tax_table(datEukr2)<-cbind(tax_table(datEukr2),labelsEuk)
tax_table(datEukr)<-cbind(tax_table(datEukr2),labelsEuk)


#sort(unique(tax_table(datEukr2)[,"Rank3"]))
#ind<-which(tax_table(datEukr2)[,"Rank3"]=="__uncultured_Eimeriidae")
#unique(tax_table(datEukr2)[ind,])
#tax_table(datEukr2)[which(rownames(tax_table(datEukr2))=="denovo64661")]

#looking where the weird centrohelid fungi comes out in a tree
ex1<-subset_taxa(datEukr2,Rank3=="__Fungi")
myTaxa<-names(taxa_sums(ex1))#[1:46]
myTaxa<-c(myTaxa,"denovo65528")
ex3<-subset_taxa(datEukr2,Rank2=="__Centrohelida")
myTaxa<-c(myTaxa,names(taxa_sums(ex3)))
ex3<-subset_taxa(datEukr2,Rank3=="__Rhizaria")
myTaxa<-c(myTaxa,names(taxa_sums(ex3)))
ex2 = prune_taxa(myTaxa, datEukr2)
plot_tree(ex2, label.tips = "taxa_names",color="Rank3")

#look at tree of abundant taxa
myTaxa<-c(names(sort(taxa_sums(datEukr2),decreasing=T)))[1:100]
ex2 = prune_taxa(myTaxa, datEukr2)
plot_tree(ex2, label.tips = "taxa_names",color="Rank3")






##### Read in ITS files, filtered with singletons removed #####
otufileITS <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/ITS/ITS_ALL_97_OTU_tablefiltsingnonchimericFunSoil2015sing.biom")
head(tax_table(otufileITS))

#Import mapping and tree file
mapITS<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/ITS_Niwot_20072015_All_MapFile.txt")

datITS<-merge_phyloseq(otufileITS,mapITS)

#rarefy and transform to relative abundance
min(sample_sums(datITS))#rarefy to 999
datITSr<-rarefy_even_depth(datITS,sample.size=min(sample_sums(datITS)),rngseed=10,replace=F) #5596 OTUs were removed because they are no longer present in any sample after random subsampling
datITSr2 = transform_sample_counts(datITSr, function(x) x/sum(x) )

head(tax_table(datITSr))
unique(tax_table(datITSr)[,"Rank2"])

#remove the k__ with substring
labelsITS<-substring(tax_table(datITSr)[,"Rank1"],4)

colnames(labelsITS)<-"labels"

#replace tax table 
tax_table(datITSr)<-cbind(tax_table(datITSr),labelsITS)

temp<-as.data.frame(sample_data(datITS))
sort(as.numeric(as.character(temp$Sample_name)))







##### Read in bacteria files, filtered with singletons removed #####
otufileBac <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Bact/16S_ALL_97_OTU_tablefiltsingnonchimerickeepbactarcfiltmitchlSoil2015readsing.biom")
head(tax_table(otufileBac))

#Import mapping and tree file
mapBac<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/515BC_Niwot_20072015_All_MapFile.txt")

treeBac<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Bact/16S_ALL_truncate_97_rep_set_filtsingS2015nonchimeras_sinaaln.tre")

datBac<-merge_phyloseq(otufileBac,mapBac,treeBac)

#rarefy and transform to relative abundance
min(sample_sums(datBac))#rarefy to 6780
datBacr<-rarefy_even_depth(datBac,sample.size=min(sample_sums(datBac)),rngseed=10,replace=F) #11966 OTUs were removed because they are no longer present in any sample after random subsampling
datBacr2 = transform_sample_counts(datBacr, function(x) x/sum(x) )

tax_table(datBacr)
unique(tax_table(datBacr)[,"Rank2"])

#remove the __ with substring
labelsBac<-substring(tax_table(datBacr)[,"Rank2"],3)

colnames(labelsBac)<-"labels"

#replace tax table 
tax_table(datBacr)<-cbind(tax_table(datBacr),labelsBac)
tax_table(datBacr2)<-cbind(tax_table(datBacr),labelsBac)

#make otu table 10:52-11:02 (takes 10 min)
datBacr2otu<-cbind(sample_data(datBacr2),t(otu_table(datBacr2)))
datBacr2otu$Sample_name<-as.numeric(as.character(datBacr2otu$Sample_name))





###### Filter data sets for network analysis ######
#Follwing Widder et al 2014 PNAS
#files: datEukr6otu, datBacr2 or datBacr2otu

#take out doubletons and singletons
datEukr7otu<-cbind(datEukr6otu[,1:31],datEukr6otu[,((which(colSums(datEukr6otu[,32:dim(datEukr6otu)[2]]>0)>2))+31)])
datBacr3otu<-cbind(datBacr2otu[,1:31],datBacr2otu[,((which(colSums(datBacr2otu[,32:dim(datBacr2otu)[2]]>0)>2))+31)])

#filter out taxa that have a summed relative abundance of <.002 (.2%)
datEukr8otu<-cbind(datEukr7otu[,1:31],datEukr7otu[,((which(colSums(datEukr7otu[,32:dim(datEukr7otu)[2]])>0.002))+31)]) #2084 otu
datBacr4otu<-cbind(datBacr3otu[,1:31],datBacr3otu[,((which(colSums(datBacr3otu[,32:dim(datBacr3otu)[2]])>0.002))+31)]) #3903 otu




#order of doing things: in qiime took out single reads (b/c likely sequencing error), in qiime filtered out unwanted taxa (chloroplasts, spiders) and samples (S.2015), in qiime took out single reads again (they are still single reads so I say they are b/c likely sequencing error), rarefied, relativized, took out doubletons and singletons, took out samples <2% summed abundance. Before I rarefied here but I think that is wrong because the doubletons/singletons/.2% otus that I removed are real. I could relativize again here, but again, the rare otus were real, they were removed only for simplification, I could have included them in the network analysis, so I won't re-relativize.























#I can't extract nematode samples OTUs because the holozoa do not have labels in the mapping file. tax_table2 is the original tax table with full line of taxononmy Eukaryota_Opisthokonta_Nucletmycea_Fungi_Dikarya_Ascomycota_Pezizomycotina_Sordariomycetes_Hypocreales
#all taxa in the file going into the network analysis (not rarefied but doesn't matter)
temp3<-rownames(tax_table(dats10)) 

#full taxonomies
head(tax_table2)
temptax<-as.character(tax_table2)
ind<-grep("Nematoda",temptax)
temp2<-tax_table2[ind,]
nematodeids<-rownames(temp2)

intersect(temp,nematodeids) #there are no nematodes present in the otufile going into network analysis
length(union(temp,nematodeids))

#checking that there are still nematodes in the original 2015 subset file, they are present in dat and dats but not dats2
temp<-rownames(otu_table(dats))
intersect(temp,nematodeids)

temp <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euk_ALL_97_OTU_filtsing.biom")
dim(tax_table(temp))
dim(otu_table(temp))

ind<-which(rownames(otu_table(temp))=="denovo244871")
tax_table(temp)[ind,] #it imported correctly and has holozoa in the label

#after my taxonomy rearrangement what did it do? that is correct
tax_table(otufile)[ind,]

#











######Grouping by kingdom#####
dats7<-transform_sample_counts(dats2, function(x) 100*x/sum(x) )
dats8<-otu_table(dats7)
dats9<-aggregate.data.frame(dats8,by=list(kingdomlabels=kingdomlabels),sum)
rownames(dats9)<-dats9$kingdomlabels
dats9$kingdomlabels<-NULL
dats9kingdom<-cbind(sample_data(dats2),t(dats9))







##### Nematode data #####
#import_biom gave an error when I tried to import the biom file I created on the server from Dorota's otu text table. I can't figure out why it gave the error beacuse the text file looks exactly like the other text files that have all euk data.
otufileN<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Nematodes/Nema_OTUtable_species.csv")
head(otufileN)
otufileN1<-otufileN[,2:(dim(otufileN)[2]-1)]
rownames(otufileN1)<-otufileN[,1]
head(otufileN1)
otufileN2 = otu_table(otufileN1, taxa_are_rows = TRUE)
taxfileN<-matrix(otufileN[,dim(otufileN)[2]])
rownames(taxfileN)<-otufileN[,1]
colnames(taxfileN)<-"Rank1"
taxfileN1 = tax_table(taxfileN)
datN = phyloseq(otufileN2,taxfileN1,map,treefile)
datN
pdN<-pd(as.matrix(as.data.frame(t(otu_table(datN)))),phy_tree(datN),include.root=TRUE) 

rownames(otu_table(otufile))[which(rownames(otu_table(otufile))%in%rownames(otufileN1))]
treefile$tip.label[which(treefile$tip.label%in%rownames(otufileN1))]
denovo120212 denovo68570
tax_table(otufile)[which(rownames(tax_table(otufile))=="denovo256418")]
#denovo256418

















#Haven't done this yet. I did this in the cooccurrencenetworkseuks file
######Read in plant data to merge with order file######
#plantcomp<-read.csv("/Users/farrer/Dropbox/Niwot Moving Uphill/Analysis/Niwot_MovingUpHill_comp2015.csv")
#head(plantcomp)
#names(plantcomp)[1]<-"Sample_name"

#Remove plant species only present in one or two plots
#dim(plantcomp)
#plantcomp2<-plantcomp[,colSums(plantcomp>0)>2]
#plantlabels<-as.data.frame(cbind(colnames(plantcomp2)[2:56],"Plant"))
#colnames(plantlabels)<-c("orders","kingdomlabels")

#Merge plants with microbes, plantcomp is everything, plantcomp2 removes doubletons/singletons
#microbplant<-merge(dats6order,plantcomp,"Sample_name",sort=F)
#microbplant2<-merge(dats6order,plantcomp2,"Sample_name",sort=F)
#head(microbplant)









