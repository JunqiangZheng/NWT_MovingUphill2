
#######Read in OTU data#######
#files:
#Euks no metazoa no fungi 97%    otufileEukS
#Euks metazoa 97%    otufileEukN
#Fungi ITS 97%
#Bacteria 97%
#Euks metazoa 99%

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


#****
#take out other samples for final reads and OTU counts
datEukSfincount<-subset_samples(datEukS,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)
datEukSfincount2<-filter_taxa(datEukSfincount, function(x) sum(x) > (0), prune=T)
sum(otu_table(datEukSfincount2))

datEukNfincount<-subset_samples(datEukN,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)
datEukNfincount2<-filter_taxa(datEukNfincount, function(x) sum(x) > (0), prune=T)
sum(otu_table(datEukNfincount2))
#****


#remove a few samples, sample 2A for nematodes (rep), sample N.2A has only 3 otus in it and they all are abundant in other samples so it doesn't affect single read count
datEukN2 <- prune_samples(sample_names(datEukN)!="N.2A.2015", datEukN)
#sample 78 for nematodes (it only had 164 reads)
datEukN3 <- prune_samples(sample_names(datEukN2)!="N.78.2015", datEukN2)
#maybe 61 which has 489 reads for the soil dataset
datEukS2 <- prune_samples(sample_names(datEukS)!="S.61.2015", datEukS)



###########for Dorota's nematode ms############
#read in file with metazoa
DEukS<-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_97_S111_OTU_tablefiltsingnonchimericbactarcplantEukSoil2015.biom")
DEukS2<-merge_phyloseq(DEukS,mapEuk)

#Filter out fungi and nematodes then calcualte Chao1 on the 18S data
which(tax_table(DEukS2)[,4]=="__Nematoda")
which(tax_table(DEukS2)[,3]=="__Fungi")
#which(tax_table(DEukS2)[,4]=="__Arthropoda")
DEukS3<-subset_taxa(DEukS2, is.na(Rank4)==T|Rank4!="__Nematoda")
DEukS4<-subset_taxa(DEukS3, is.na(Rank3)==T|Rank3!="__Fungi")

#filter singletons
DEukS5<-prune_taxa(taxa_sums(DEukS4) > 1, DEukS4)

#filter plot with only 545 reads
DEukS6 <- prune_samples(sample_names(DEukS5)!="S.61.2015", DEukS5)

#rarefy
sort(sample_sums(DEukS6))#rarefy to 939
DEukS7<-rarefy_even_depth(DEukS6,sample.size=min(sample_sums(DEukS6)),rngseed=10,replace=F) #

#calculate Chao1
DChao1<-estimate_richness(DEukS7, measures=c("Chao1", "Shannon","Observed"))
DChao1a<-DChao1[order(rownames(DChao1)),]
write.csv(DChao1a,"/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/18SChao1nofungiornematodesforDorotaNEW.csv",row.names=T)
#####################


#rarefy and transform to relative abundance
min(sample_sums(datEukS2))#rarefy to 808
min(sample_sums(datEukN3))#rarefy to 692
datEukS3<-rarefy_even_depth(datEukS2,sample.size=min(sample_sums(datEukS2)),rngseed=10,replace=F) #1895 OTUs were removed because they are no longer present in any sample after random subsampling, 
datEukN4<-rarefy_even_depth(datEukN3,sample.size=min(sample_sums(datEukN3)),rngseed=10,replace=F) #2601 OTUs were removed because they are no longer present in any sample after random subsampling, 
datEukS4 = transform_sample_counts(datEukS3, function(x) x/sum(x) )
datEukN5 = transform_sample_counts(datEukN4, function(x) x/sum(x) )



#make a label column with kingdom/phylum level labels

#Making groups for labeling graphs
#Amoebozoa (kingdom)
#photosynthetic Alveolata (kingdom) (phylum Dinoflagellata: mostly photosynthetic; nonphotosynthetic Protoperidinium, both SL163A10 and Pfiesteria (can if it eats an alga), unknown D244), 
#nonphotosynthetic Alveolata (phyla Ciliophora(predators), protalveolata, apicomplexa (parasites)), BOLA553 is a near apicomplexan, not sure what SCM37C52 is so I put it here
#Archaeplastida (major_clade), combine the kingdoms Chloroplastida, Rhodophyceae, Haptophyta, 
#Rhizaria (kingdom: unicellular amoeboid euks), 
#Holozoa(kingdom: all animals, not fungi), includes metazoa (animals)
#photosynthetic Excavata major clade (was Discoba kingdom) (class Euglenida: mostly photosynthetic), 
#nonphotosnthetic Excavata (was Discoba) (in discoba, phylum Heterolobosea: parasites, free living, symbiotic, amoeba-like, and phylum Jakoba heterotrophic), superphylum metamonada, doesn't have kingdom, parasites not Ps.
#photosynthetic Cryptophyceae is a kingdom, they have one or two chloroplasts, except for Chilomonas, which has leucoplasts (not pigmented not for Ps) and Goniomonas (formerly Cyathomonas) which lacks plastids entirely.P1-31 (i can't find much info on this, some sites called them picoplankton and most of the cryptophyceae are photosynthetic). Cryptomonadales can't find much about it but assume Ps. __SA1-3C06 is at rank2 but notes in silva say "cryptophyte"
#nonphotosynthetic cryptophyceae - Chilomonas and Goniomonas (formerly Cyathomonas)
#Haptophyta kingdom, I think all are Ps
#Centrohelida - encyclopedia of life said some species have photosynthetic symbionts but I can't find any info on Ps taxa.Heterophryidae, Mb-5C (can't find any info on, I'll put it with nonPs), M1-18D08(can't find info),Acanthocystidae(can't find info),Pterocystis(can't find info), so I will leave them all in a group called "centrohelida"
#NonPs Stramenopiles: MAST-3, MAST-12, MAST-7, Labyrinthulomycetes,Bicosoecida,Peronosporomycetes,Hyphochytriales,__Incertae_Sedis(were all __Pirsonia which is a parasite of algae)
#Ps Stramenopiles: Pelagophycea, Diatomea, Eustigmatales, Xanthophyceae, Phaeothamniophyceae, Chrysophyceae,Ochromonas,CCI40 (not much info, some places call it is chrysophyte),Synurales, Raphidophyceae, NA1-2A5 (can't find any info so putting it with Ps), TKR07M.92

#heterotrophic eukaryota (things without a kingdom): Picozoa (phylum nonphotosynthetic unicellular eukaryotes, only 1 OTU); RT5iin14; Apusomonadidae (nonphotosynthetic protist, unknown kingdom); Kathablepharidae (nonPs protist, unkonwn kingdo);Telonema (nonPs protist genus);Breviatea (amoeba like);DH147-EKD23;__Ancyromonadida;GoC2-B10 (no info);Zeuk77 (no info)

#datEukS4
#datEukN5

head(tax_table(datEukS4))
labelsEukS<-tax_table(datEukS4)[,"Rank3"]#

ind<-which(tax_table(datEukS4)[,"Rank2"]=="__Amoebozoa")
labelsEukS[ind]<-"Amoebozoa"
ind<-which(tax_table(datEukS4)[,"Rank4"]=="__Dinoflagellata")#__Protoperidinium is included but I then change it below becuase some things have an NA for rank8
labelsEukS[ind]<-"Photosynthetic_Alveolata"
ind<-which(tax_table(datEukS4)[,"Rank8"]=="__Protoperidinium")
labelsEukS[ind]<-"Nonphotosynthetic_Alveolata"
ind<-which(tax_table(datEukS4)[,"Rank4"]!="__Dinoflagellata"&tax_table(datEukS4)[,"Rank3"]=="__Alveolata")
labelsEukS[ind]<-"Nonphotosynthetic_Alveolata"
ind<-which(tax_table(datEukS4)[,"Rank2"]=="__Archaeplastida")
labelsEukS[ind]<-"Archaeplastida"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Rhizaria")
labelsEukS[ind]<-"Rhizaria"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Holozoa")
labelsEukS[ind]<-"Holozoa"
ind<-which(tax_table(datEukS4)[,"Rank6"]=="__Euglenida")
labelsEukS[ind]<-"Photosynthetic_Excavata"
ind<-which(tax_table(datEukS4)[,"Rank6"]!="__Euglenida"&tax_table(datEukS4)[,"Rank2"]=="__Excavata")
labelsEukS[ind]<-"Nonphotosynthetic_Excavata"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Goniomonas")
labelsEukS[ind]<-"Nonphotosynthetic_Chryptophyceae"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__P1-31"|tax_table(datEukS4)[,"Rank3"]=="__Cryptomonadales"|tax_table(datEukS4)[,"Rank2"]=="__SA1-3C06")
labelsEukS[ind]<-"Photosynthetic_Chryptophyceae"
ind<-which(tax_table(datEukS4)[,"Rank2"]=="__Haptophyta")
labelsEukS[ind]<-"Haptophyta"
ind<-which(tax_table(datEukS4)[,"Rank2"]=="__Centrohelida"&tax_table(datEukS4)[,"Rank3"]!="__uncultured_Sarcosomataceae")
labelsEukS[ind]<-"Centrohelida"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__uncultured_Sarcosomataceae")#i'm putting the wierd centrohelid/fungus with centrohilid (since I supposedly took out al fungi)
labelsEukS[ind]<-"Centrohelida"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Stramenopiles"&tax_table(datEukS4)[,"Rank4"]%in%c("__MAST-3","__MAST-12","__MAST-7","__Labyrinthulomycetes","__Bicosoecida","__Peronosporomycetes","__Hyphochytriales","__Incertae_Sedis"))
labelsEukS[ind]<-"Nonphotosynthetic_Stramenopiles"
ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Stramenopiles"&tax_table(datEukS4)[,"Rank4"]%in%c("__Pelagophycea","__Diatomea","__Eustigmatales","__Xanthophyceae","__Phaeothamniophyceae","__Chrysophyceae","__Ochromonas","__CCI40","__Synurales","__Raphidophyceae","__NA1-2A5","__TKR07M.92","__Pelagophyceae"))
labelsEukS[ind]<-"Photosynthetic_Stramenopiles"

ind<-which(tax_table(datEukS4)[,"Rank3"]=="__RT5iin14"|tax_table(datEukS4)[,"Rank2"]=="__Picozoa"|tax_table(datEukS4)[,"Rank2"]=="__RT5iin25"|tax_table(datEukS4)[,"Rank3"]=="__Apusomonadidae"|tax_table(datEukS4)[,"Rank2"]=="__Kathablepharidae"|tax_table(datEukS4)[,"Rank3"]=="__Telonema"|tax_table(datEukS4)[,"Rank3"]=="__Breviatea"|tax_table(datEukS4)[,"Rank2"]=="__DH147-EKD23"|tax_table(datEukS4)[,"Rank3"]=="__Ancyromonadida"|tax_table(datEukS4)[,"Rank2"]=="__GoC2-B10"|tax_table(datEukS4)[,"Rank2"]=="__Zeuk77")
labelsEukS[ind]<-"Heterotrophic_Eukarya"

#unique(labelsEukS)
#sort(unique(tax_table(datEukS4)[,"Rank3"]))
#ind<-which(tax_table(datEukS4)[,"Rank3"]=="__Alveolata")
#unique(tax_table(datEukS4)[ind,"Rank4"])
#tax_table(datEukS4)[which(rownames(tax_table(datEukS4))=="denovo64661")]

#I could separate excavata into kingdoms/super phyla and also separate archaeplastida into its 3 kingdoms
unique(labelsEukS)
colnames(labelsEukS)<-"labels"

#replace tax table for datEukS4 (relative abun) and datEukS3 (not relativized)
tax_table(datEukS4)<-cbind(tax_table(datEukS4),labelsEukS)
tax_table(datEukS3)<-cbind(tax_table(datEukS3),labelsEukS)

#looking where the weird centrohelid fungi comes out in a tree
#ex1<-subset_taxa(datEukr2,Rank3=="__Fungi")
#myTaxa<-names(taxa_sums(ex1))#[1:46]
# myTaxa<-c(myTaxa,"denovo65528")
# ex3<-subset_taxa(datEukr2,Rank2=="__Centrohelida")
# myTaxa<-c(myTaxa,names(taxa_sums(ex3)))
# ex3<-subset_taxa(datEukr2,Rank3=="__Rhizaria")
# myTaxa<-c(myTaxa,names(taxa_sums(ex3)))
# ex2 = prune_taxa(myTaxa, datEukr2)
# plot_tree(ex2, label.tips = "taxa_names",color="Rank3")

#look at tree of abundant taxa
#myTaxa<-c(names(sort(taxa_sums(datEukr2),decreasing=T)))[1:100]
#ex2 = prune_taxa(myTaxa, datEukr2)
#plot_tree(ex2, label.tips = "taxa_names",color="Rank3")


#labels for datEukN5 (relativized) and datEukN4 (not relativized), use Rank4
head(tax_table(datEukN5))
labelsEukN<-substring(tax_table(datEukN5)[,"Rank4"],3)
colnames(labelsEukN)<-"labels"

#replace tax table 
tax_table(datEukN5)<-cbind(tax_table(datEukN5),labelsEukN)
tax_table(datEukN4)<-cbind(tax_table(datEukN4),labelsEukN)







##### Read in ITS files, filtered with singletons removed #####
otufileITS <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/ITS/ITS_ALL_97_OTU_tablefiltsingnonchimericFunSoil2015sing.biom")
head(tax_table(otufileITS))

#Import mapping and tree file
mapITS<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/ITS_Niwot_20072015_All_MapFilenewlomehi.txt")

datITS<-merge_phyloseq(otufileITS,mapITS)

#****
#take out other samples for final reads and OTU counts
datITSfincount<-subset_samples(datITS,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)
datITSfincount2<-filter_taxa(datITSfincount, function(x) sum(x) > (0), prune=T)
sum(otu_table(datITSfincount2))
#****

#rarefy and transform to relative abundance
min(sample_sums(datITS))#rarefy to 999
datITS2<-rarefy_even_depth(datITS,sample.size=min(sample_sums(datITS)),rngseed=10,replace=F) #5596 OTUs were removed because they are no longer present in any sample after random subsampling
datITS3 = transform_sample_counts(datITS2, function(x) x/sum(x) )

head(tax_table(datITS3))
unique(tax_table(datITS3)[,"Rank2"])

#remove the p__ with substring
labelsITS<-substring(tax_table(datITS3)[,"Rank2"],4)

colnames(labelsITS)<-"labels"

#replace tax table, datITS3 (relativized) datITS2 (not relativized)
tax_table(datITS3)<-cbind(tax_table(datITS3),labelsITS)
tax_table(datITS2)<-cbind(tax_table(datITS2),labelsITS)








##### Read in bacteria files, filtered with singletons removed #####
otufileBac <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Bact/16S_ALL_97_OTU_tablefiltsingnonchimerickeepbactarcfiltmitchlSoil2015readsing.biom")
head(tax_table(otufileBac))

#Import mapping and tree file
mapBac<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/515BC_Niwot_20072015_All_MapFilenewlomehi.txt")

treeBac<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Bact/16S_ALL_truncate_97_rep_set_filtsingS2015nonchimeras_sinaaln.tre")

datBac<-merge_phyloseq(otufileBac,mapBac,treeBac)

#****
#take out other samples for final reads and OTU counts
datBacfincount<-subset_samples(datBac,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)
datBacfincount2<-filter_taxa(datBacfincount, function(x) sum(x) > (0), prune=T)
sum(otu_table(datBacfincount2))
#****

#rarefy and transform to relative abundance
min(sample_sums(datBac))#rarefy to 6780
datBac2<-rarefy_even_depth(datBac,sample.size=min(sample_sums(datBac)),rngseed=10,replace=F) #11966 OTUs were removed because they are no longer present in any sample after random subsampling
datBac3 = transform_sample_counts(datBac2, function(x) x/sum(x) )

tax_table(datBac3)
unique(tax_table(datBac3)[,"Rank2"])

#remove the __ with substring
labelsBac<-substring(tax_table(datBac3)[,"Rank2"],3)

#make the class the label for anything in the Phylum Chloroflexi, Only the class Chloroflexi are photosynthetic, the other classes (Ktedontobacteria) are not photosynthetic. There many OTU listed as phylum Chloroflexi but that do not have a Class. The majority of the taxa in Chloroflexi that do have a class are Ktedonobacteria. Thus I will make the default color for heterotrophs (Purple), and just change the Chloroflexales to photosynthetic (green)
ind<-which(labelsBac=="Chloroflexi")
temp<-data.frame(tax_table(datBac3))
length(which(temp$Rank2=="__Chloroflexi"))
length(which(temp$Rank2=="__Chloroflexi"&temp$Rank3=="__Chloroflexi"))#none are known to be in class chloroflexi
length(which(temp$Rank3=="__Chloroflexales"|temp$Rank4=="__Chloroflexaceae")) #81 are Chloroflexales which I imagine are in the class Chloroflexi
length(which(temp$Rank2=="__Chloroflexi"&temp$Rank4=="__Chloroflexaceae")) #61 are Chloroflexaceae, these are included in the 81 above
length(which(temp$Rank2=="__Chloroflexi"&temp$Rank3=="__Ktedonobacteria"))#1497 are Ktedonobacteria
temp<-subset(temp,Rank2=="__Chloroflexi")
sort(unique(temp$Rank3))
sort(unique(temp$Rank4))
#to label the chloroflexales as phototrophs
temp<-data.frame(tax_table(datBac3))
ind<-which(temp$Rank3=="__Chloroflexales"|temp$Rank4=="__Chloroflexaceae")
labelsBacchloro<-labelsBac
labelsBacchloro[ind]<-"Chloroflexiphoto"
#to label the ktedonobacteria as heterotrophs
# temp<-data.frame(tax_table(datBac3))
# ind<-which(temp$Rank3=="__Ktedonobacteria")
# labelsBacchloro<-labelsBac
# labelsBacchloro[ind]<-"Chloroflexihetero"

colnames(labelsBac)<-"labels"
colnames(labelsBacchloro)<-"labels"

#replace tax table, datBac3 (relativized) datBac2 (not relativized)
tax_table(datBac2)<-cbind(tax_table(datBac2),labelsBac)
tax_table(datBac3)<-cbind(tax_table(datBac3),labelsBac)





###### Take out samples that aren't shared across all datasets ######
#For euksS sample 81 did not amplify and 61 had low # reads. for euksN sample 33 and 56 did not have enough soil so were not done, and 78 had low # reads. for ITS 126 didnt amplify. for bacteria, samples 5,34,126 did not amplify. should have 90 samples left
#Files: datBac2, datBac3, datEukN4, datEukN5, datEukS3, datEukS4, datITS2, datITS3

#notes: it looks like both prune_samples and subset_samples do NOT remove taxa that have an abundance of 0 after the samples are removed.
#"%w/o%" <- function(x, y) x[!x %in% y] #that function is not exactly what i wanted, it subsets the list rather than returning indices. this is how I would do it which(!test%in%notlist) except that I can't seem to access the sample_data in the prune_samples code
#notlist<-c(81,61,33,56,78,126,5,34)

datBac2f <- subset_samples(datBac2,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)
datBac3f <- subset_samples(datBac3,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)

datEukN4f <- subset_samples(datEukN4,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)
datEukN5f <- subset_samples(datEukN5,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)

datEukS3f <- subset_samples(datEukS3,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)
datEukS4f <- subset_samples(datEukS4,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)

datITS2f <- subset_samples(datITS2,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)
datITS3f <- subset_samples(datITS3,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)


###I should have taken out taxa that are now zero across all samples due to deleting above the one sample that has the taxon. filter_taxa(datITS3f, function(x) sum(x) > (0), prune=T)


#make otu tables (bact takes 10 min)
datBacr3fotu<-cbind(sample_data(datBac3f),t(otu_table(datBac3f)))
datBacr3fotu$Sample_name<-as.numeric(as.character(datBacr3fotu$Sample_name))

datEukN5fotu<-cbind(sample_data(datEukN5f),t(otu_table(datEukN5f)))
datEukN5fotu$Sample_name<-as.numeric(as.character(datEukN5fotu$Sample_name))

datEukS4fotu<-cbind(sample_data(datEukS4f),t(otu_table(datEukS4f)))
datEukS4fotu$Sample_name<-as.numeric(as.character(datEukS4fotu$Sample_name))

datITS3fotu<-cbind(sample_data(datITS3f),t(otu_table(datITS3f)))
datITS3fotu$Sample_name<-as.numeric(as.character(datITS3fotu$Sample_name))




##### Phylogenetic diversity ######
pdBac<-pd(as.matrix(datBacr3fotu[,-c(1:31)]),phy_tree(datBac3f),include.root=TRUE) #takes a few hours, started 3:45pm, finished ~midnight
pdEukN<-pd(as.matrix(datEukN5fotu[,-c(1:31)]),phy_tree(datEukN5f),include.root=TRUE) #took 3 minutes
pdEukS<-pd(as.matrix(datEukS4fotu[,-c(1:31)]),phy_tree(datEukS4f),include.root=TRUE) #took 5 minutes
richITS<-as.data.frame(rowSums(datITS3fotu[,-c(1:31)]>0))
pdMet<-pd(as.matrix(dat99Met2fotu[,-c(1:31)]),phy_tree(dat99Met2f),include.root=TRUE) #run code from metazoa section below first
pdNem<-pd(as.matrix(dat99Met2fNemotu[,-c(1:31)]),phy_tree(dat99Met2fNem),include.root=TRUE) # I need to use root=T, if I use root=F it cannot calculate the diversity in a sample with only one taxon



###### Grouping by kingdom/phylum #####

datBac3fk<-aggregate.data.frame(otu_table(datBac3f),by=list(labels=tax_table(datBac3f)[,"labels"]),sum)
rownames(datBac3fk)<-datBac3fk$labels
datBac3fk$labels<-NULL
datBac3fk2<-cbind(sample_data(datBac3f),t(datBac3fk))
head(datBac3fk2)

datEukN5fk<-aggregate.data.frame(otu_table(datEukN5f),by=list(labels=tax_table(datEukN5f)[,"labels"]),sum)
rownames(datEukN5fk)<-datEukN5fk$labels
datEukN5fk$labels<-NULL
datEukN5fk2<-cbind(sample_data(datEukN5f),t(datEukN5fk))
head(datEukN5fk2)

datEukS4fk<-aggregate.data.frame(otu_table(datEukS4f),by=list(labels=tax_table(datEukS4f)[,"labels"]),sum)
rownames(datEukS4fk)<-datEukS4fk$labels
datEukS4fk$labels<-NULL
datEukS4fk2<-cbind(sample_data(datEukS4f),t(datEukS4fk))
head(datEukS4fk2)

datITS3fk<-aggregate.data.frame(otu_table(datITS3f),by=list(labels=tax_table(datITS3f)[,"labels"]),sum)
rownames(datITS3fk)<-datITS3fk$labels
datITS3fk$labels<-NULL
datITS3fk2<-cbind(sample_data(datITS3f),t(datITS3fk))
head(datITS3fk2)


dat99Met2fNem2<-transform_sample_counts(dat99Met2fNem, function(x) x/sum(x) )
temp<-otu_table(dat99Met2fNem2)
temp2<-temp[order(rownames(temp)),]
labels99Nemr2<-labels99Nemr[order(rownames(labels99Nemr)),]
dat99Met2fNem2k<-aggregate.data.frame(temp2,by=list(labels=labels99Nemr2$labels),sum)
rownames(dat99Met2fNem2k)<-dat99Met2fNem2k$labels
dat99Met2fNem2k$labels<-NULL
dat99Met2fNem2k2<-cbind(sample_data(dat99Met2fNem2),t(dat99Met2fNem2k))
colnames(dat99Met2fNem2k2)[32:38]<-c("Animal_feeder","Animal_parasite","Bacterial_feeder","Fungal_feeder","Omnivore","Plant_parasite","Root_associate")




###### Filter data sets for network analysis ######
#Follwing Widder et al 2014 PNAS
#files: datBacr3fotu datEukN5fotu datEukS4fotu datITS3fotu

#take out doubletons and singletons
datBacr3fotu2<-cbind(datBacr3fotu[,1:31],datBacr3fotu[,((which(colSums(datBacr3fotu[,32:dim(datBacr3fotu)[2]]>0)>2))+31)])

datEukN5fotu2<-cbind(datEukN5fotu[,1:31],datEukN5fotu[,((which(colSums(datEukN5fotu[,32:dim(datEukN5fotu)[2]]>0)>2))+31)])

datEukS4fotu2<-cbind(datEukS4fotu[,1:31],datEukS4fotu[,((which(colSums(datEukS4fotu[,32:dim(datEukS4fotu)[2]]>0)>2))+31)])

datITS3fotu2<-cbind(datITS3fotu[,1:31],datITS3fotu[,((which(colSums(datITS3fotu[,32:dim(datITS3fotu)[2]]>0)>2))+31)])

#filter out taxa that have a summed relative abundance of <.002 (.2%)
datBacr3fotu3<-cbind(datBacr3fotu2[,1:31],datBacr3fotu2[,((which(colSums(datBacr3fotu2[,32:dim(datBacr3fotu2)[2]])>0.002))+31)]) #3772 otu

datEukN5fotu3<-cbind(datEukN5fotu2[,1:31],datEukN5fotu2[,((which(colSums(datEukN5fotu2[,32:dim(datEukN5fotu2)[2]])>0.002))+31)]) #233 otu

datEukS4fotu3<-cbind(datEukS4fotu2[,1:31],datEukS4fotu2[,((which(colSums(datEukS4fotu2[,32:dim(datEukS4fotu2)[2]])>0.002))+31)]) #1498 otu

datITS3fotu3<-cbind(datITS3fotu2[,1:31],datITS3fotu2[,((which(colSums(datITS3fotu2[,32:dim(datITS3fotu2)[2]])>0.002))+31)]) #1125 otu



#order of doing things: in qiime took out single reads (b/c likely sequencing error), in qiime filtered out unwanted taxa (chloroplasts, spiders) and samples (S.2015), in qiime took out single reads again (they are still single reads so I say they are b/c likely sequencing error), rarefied, relativized, took out doubletons and singletons, took out samples <2% summed abundance. Before I rarefied here but I think that is wrong because the doubletons/singletons/.2% otus that I removed are real. I could relativize again here, but again, the rare otus were real, they were removed only for simplification, I could have included them in the network analysis, so I won't re-relativize.






###### Make labelfile for network analysis ######

#test if there is denovo name overlap, yes a ton, so relabel the denovos here and in the labels files
namesBac<-names(datBacr3fotu3[,-c(1:31)])
namesEukN<-names(datEukN5fotu3[,-c(1:31)])
namesEukS<-names(datEukS4fotu3[,-c(1:31)])
namesITS<-names(datITS3fotu3[,-c(1:31)])
length(namesBac)
length(namesITS)
length(union(namesBac,namesITS))
intersect(namesBac,namesITS)

namesBac2 <- sub("^", "b", namesBac)
namesEukN2 <- sub("^", "n", namesEukN)
namesEukS2 <- sub("^", "s", namesEukS)
namesITS2 <- sub("^", "i", namesITS)

names(datBacr3fotu3)[-c(1:31)]<-namesBac2
names(datEukN5fotu3)[-c(1:31)]<-namesEukN2
names(datEukS4fotu3)[-c(1:31)]<-namesEukS2
names(datITS3fotu3)[-c(1:31)]<-namesITS2

#I have tow labels files for bacteria, the labelsBac is by phylum, the labelsBac2is for networks where some of the chloroflexi are heterotrophs and some are phototrophs
labelsBac2<-labelsBac
rownames(labelsBac2)<-sub("^", "b",rownames(labelsBac2))
labelsBacchloro2<-labelsBacchloro
rownames(labelsBacchloro2)<-sub("^", "b",rownames(labelsBacchloro2))
labelsEukN2<-labelsEukN
rownames(labelsEukN2)<-sub("^", "n",rownames(labelsEukN2))
labelsEukS2<-labelsEukS
rownames(labelsEukS2)<-sub("^", "s",rownames(labelsEukS2))
labelsITS2<-labelsITS
rownames(labelsITS2)<-sub("^", "i",rownames(labelsITS2))










##### Read in euk 99% metazoa files, counts, from Dorota's cleanup #####
#I read in all metazoa, but then below I subset out only the nematodes
otufile99Met <-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_AllMeatazoa_99_S11_OTU.csv")
#it is not importing the biom file into phyloseq, I had this import problem before, when I tried to create a biom from txt file. So I will have to use the txt/csv file and then try to merge that
#library("biom")
#biom = read_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_AllMeatazoa_99_S11_OTU.biom") #this works but still cannot be merged with mapping files in phyloseq

otufile99Met2<-as.matrix(otufile99Met[,2:97])
rownames(otufile99Met2)<-otufile99Met[,1]
taxtable99Met<-as.character(otufile99Met[,98])

mylist<-strsplit(taxtable99Met, "; ", TRUE)
mymatrix<-matrix(NA,nrow=268,ncol=8)
for(i in 1:268){
  mymatrix[i,1:length(mylist[[i]])]<-mylist[[i]]
}
colnames(mymatrix)<-c("Rank1","Rank2","Rank3","Rank4","Rank5","Rank6","Rank7","Rank8")
rownames(mymatrix)<-otufile99Met[,1]

taxtable99Met2<-mymatrix

otufile99Met3 = otu_table(otufile99Met2, taxa_are_rows = TRUE)
taxtable99Met3 = tax_table(taxtable99Met2)


#Mapping ifle is already imported 
mapEuk

#Read in tree file
tree99Met<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_AllMetazoa_truncate_99_sinaal.tre")

dat99Met<-merge_phyloseq(otufile99Met3,taxtable99Met3,mapEuk,tree99Met)

#remove craniata, there are no NAs so subset is fine, unique(tax_table(dat99Met2)[,"Rank4"])
dat99Met2 = subset_taxa(dat99Met, Rank4 != "__Craniata")


#****
#take out other samples for final reads and OTU counts
dat99Metfincount<-subset_samples(dat99Met2,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)
dat99Metfincount2<-filter_taxa(dat99Metfincount, function(x) sum(x) > (0), prune=T)
sum(otu_table(dat99Metfincount2))
#****


labels99Met<-substring(tax_table(dat99Met2)[,"Rank4"],3)

colnames(labels99Met)<-"labels"

#replace tax table in dat99Met2
tax_table(dat99Met2)<-cbind(tax_table(dat99Met2),labels99Met)

#final dataset for all metazoa
dat99Met2f <- subset_samples(dat99Met2,Sample_name!=81&Sample_name!=61&Sample_name!=33&Sample_name!=56&Sample_name!=78&Sample_name!=126&Sample_name!=5&Sample_name!=34)

#final dataset for nematodes only
dat99Met2fNem = subset_taxa(dat99Met2f, Rank4 == "__Nematoda")
sum(otu_table(dat99Met2fNem)) #final count for nematodes only

#make otu tables
dat99Met2fotu<-cbind(sample_data(dat99Met2f),t(otu_table(dat99Met2f)))
dat99Met2fotu$Sample_name<-as.numeric(as.character(dat99Met2fotu$Sample_name))

dat99Met2fNemotu<-cbind(sample_data(dat99Met2fNem),t(otu_table(dat99Met2fNem)))
dat99Met2fNemotu$Sample_name<-as.numeric(as.character(dat99Met2fNemotu$Sample_name))


###### Filter for network analysis ######
#take out doubletons and singletons
dat99Met2fotu2<-cbind(dat99Met2fotu[,1:31],dat99Met2fotu[,((which(colSums(dat99Met2fotu[,32:dim(dat99Met2fotu)[2]]>0)>2))+31)])

#filter out taxa that have a summed relative abundance of <.002 (.2%)
#transform to relative abundance
dat99Met2fotu3 <- cbind(dat99Met2fotu2[,1:31],dat99Met2fotu2[32:129]/rowSums(dat99Met2fotu2[32:129]))
min(colSums(dat99Met2fotu3[,32:129],na.rm=T))
#all of the taxa have at least a rel abun of .2% (the smallest is 0.002088). so this step is not needed

dat99Met2fotu2 #98 otu of all metazoa


#Make label files and rename denovos in otutables
names99Met<-names(dat99Met2fotu2[,-c(1:31)])

names99Met2 <- sub("^", "m", names99Met)

names(dat99Met2fotu2)[-c(1:31)]<-names99Met2

labels99Met2<-labels99Met
rownames(labels99Met2)<-sub("^", "m",rownames(labels99Met2))

head(labels99Met2)

#for manuscript, because I did the networks with all metazoa and then filtered out everything except nematodes. so how many nematodes ended up going into that analysis
temp<-merge(data.frame(otu=rownames(labels99Met2),labels99Met2),data.frame(otu=colnames(dat99Met2fotu2[,-c(1:31)])))
sum(temp$labels=="Nematoda") #63 nematode have a frequency >=3


#Labels for trophic groups for nematodes 99% only
trophicgroup<-read.csv("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_All_Nema_99_S111functionalgroup.csv")
head(trophicgroup)

#for network figures
labels99Nem<-data.frame(labels=trophicgroup[,3])
rownames(labels99Nem)<-trophicgroup[,1]
rownames(labels99Nem)<-sub("^", "m",rownames(labels99Nem))

#for relative abundance figures
labels99Nemr<-data.frame(labels=trophicgroup[,3])
rownames(labels99Nemr)<-trophicgroup[,1]
labels99Nemr$otu<-trophicgroup[,1]

#adding m labels to trophic group data.frame
trophicgroup$otu2<-rownames(labels99Nem)

