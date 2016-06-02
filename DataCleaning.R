
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

#rarefy and transform to relative abundance
min(sample_sums(datITS))#rarefy to 999
datITS2<-rarefy_even_depth(datITS,sample.size=min(sample_sums(datITS)),rngseed=10,replace=F) #5596 OTUs were removed because they are no longer present in any sample after random subsampling
datITS3 = transform_sample_counts(datITS2, function(x) x/sum(x) )

head(tax_table(datITS3))
unique(tax_table(datITS3)[,"Rank2"])

#remove the k__ with substring
labelsITS<-substring(tax_table(datITS3)[,"Rank1"],4)

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

#rarefy and transform to relative abundance
min(sample_sums(datBac))#rarefy to 6780
datBac2<-rarefy_even_depth(datBac,sample.size=min(sample_sums(datBac)),rngseed=10,replace=F) #11966 OTUs were removed because they are no longer present in any sample after random subsampling
datBac3 = transform_sample_counts(datBac2, function(x) x/sum(x) )

tax_table(datBac3)
unique(tax_table(datBac3)[,"Rank2"])

#remove the __ with substring
labelsBac<-substring(tax_table(datBac3)[,"Rank2"],3)

colnames(labelsBac)<-"labels"

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



#make otu tables (bact takes 10 min)
datBacr3fotu<-cbind(sample_data(datBac3f),t(otu_table(datBac3f)))
datBacr3fotu$Sample_name<-as.numeric(as.character(datBacr3fotu$Sample_name))

datEukN5fotu<-cbind(sample_data(datEukN5f),t(otu_table(datEukN5f)))
datEukN5fotu$Sample_name<-as.numeric(as.character(datEukN5fotu$Sample_name))

datEukS4fotu<-cbind(sample_data(datEukS4f),t(otu_table(datEukS4f)))
datEukS4fotu$Sample_name<-as.numeric(as.character(datEukS4fotu$Sample_name))

datITS3fotu<-cbind(sample_data(datITS3f),t(otu_table(datITS3f)))
datITS3fotu$Sample_name<-as.numeric(as.character(datITS3fotu$Sample_name))





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











###have not done below yet




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









