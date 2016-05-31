
#######Read in OTU data#######

##### Read in euk files, filtered with singletons removed #####
otufileEuk <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_97_S111_OTU_tablefiltsingnonchimericbactarcplantMEukSoil2015sing.biom")
#otufileEuknoFungi <-import_biom("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_97_S111_OTU_tablefiltsingnonchimericbactarcplantEukSoil2015singnoFungi.biom")
head(tax_table(otufileEuk))

#Import mapping and tree file
mapEuk<-import_qiime_sample_data("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MappingFiles/EukBr_Niwot_20072015_All_MapFile.txt")

treeEuk<-read_tree("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/Euks/Euk_ALL_truncate_97_rep_set_filtsingS2015nonchimeras_sinaaln.tre")

datEuk<-merge_phyloseq(otufileEuk,mapEuk,treeEuk)

#what other macrobes to remove from Euk dataset?
#To remove
Rank4=__Acanthocephala
Rank4=__Annelida
Rank4=__Bryozoa
Rank4=__Chaetognatha
Rank4=__Cnidaria
Rank4=__Craniata
Rank4=__Ctenophora
Rank4=__Echinodermata
Rank4=__Mollusca
Rank4=__Nemertea
Rank4=__Platyhelminthes
Rank4=__Porifera
Rank6=__Merostomata (__Arthropoda; __Chelicerata; )
Rank6= __Malacostraca(__Arthropoda; __Crustacea;)
Rank5=__Hexapoda (__Arthropoda) (all of hexapods, most are insects and collembola but even proturans at .5-2mm)
Rank5=__Myriapoda (__Arthropoda)

Arthopods, to remove (Rank 7): 
Rank7=__Eusimonia (spider); 
Rank7=__Microcaeculus (mites but the pics look larger than 1mm, maybe 2mm);
Rank7=__Neogovea (mite but looks about 2mm); 
Rank7=__Penthaleus (mites but large); 
Rank7=__Pimoa (spider);
Rank7=__Rhagidiidae (mite but look large not much info);
Rank7=__Sisicottus; (spider);
Rank7=__uncultured_eukaryote (I guess delete these)



To Keep:
Arthopods to keep:
rank7=__Alicorhagia (mites, prob less than 1mm); __Alicorhagidia (mites prob less than 1mm), __Bimichaelia (mites, prob small), __Cercoleipus (mites prob small); __Cerophagus (mites prob small); __Cosmolaelaps; (mite, prob small); __Dactylothrombium (1-1.2mm), __Derocheilocaris_typicus (.5 mm really a crustacean but it is under arachnida here); __Eupelops (~1mm);__Eupodidae (mites, most likely <1mm);__Eutegaeus (mite, prob small); __Gamasiphis (mite ~1mm); __Gehypochthonius (mite, ~1mm); __Gozmanyina (mite, cant find info about size);__Hydrozetes (1mm); __Nanhermannia (mites 1mm); __Nesiacarus (mite, no info);__Perlohmannia (about 1mm maybe a tad larger);__Rhizoglyphus (mites small); __Rostrozetes (mites small); __Tectocepheus (mite small); __Tetranychus (spider mite, .4mm but pics look larger); __Trhypochthonius (about 1mm); __Xenillus (mite about 1mm); __Zercon (mite, not much info possibly small)
Metazoans to keep: 
Rank4=nematoda, rotifera, tardigrada (.5mm),  __Gastrotricha(small worms .3mm), __Kinorhyncha (mud dragons max 1mm), __Loricifera (.1-1mm), __Placozoa (.5 or 1mm)
rank6=__Maxillopoda(__Arthropoda; __Crustacea; )
Rank6=__Ostracoda(__Arthropoda; __Crustacea)

#if I want to do it here, it actually is easier in qiime since you don't have to specify rank
#datEuk<-subset_taxa(datEuk,domain=="Eukaryota"&class!="Embryophyta")



#rarefy and transform to relative abundance
min(sample_sums(datEuk))#rarefy to 1329, 1205
datEukr<-rarefy_even_depth(datEuk,sample.size=min(sample_sums(datEuk)),rngseed=10,replace=F) #4163, 3606 OTUs were removed because they are no longer present in any sample after random subsampling, 
datEukr2 = transform_sample_counts(datEukr, function(x) x/sum(x) )
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

head(tax_table(datEukr2))
labelsEuk<-tax_table(datEukr2)[,"Rank3"]#

ind<-which(tax_table(datEukr2)[,"Rank3"]=="__Fungi")
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






###### Replace fungi from 18S with ITS ######

datITSr2
datEukr2

head(sample_data(datEukr2))

#figure out which samples are shared across Euk and ITS datasets. Euk doesnt have 81, ITS doesn't have 126 (no reads amplified at all)
temp<-as.data.frame(as.matrix(sample_data(datEukr2)[,"Sample_name"]))
sort(as.numeric(temp$Sample_name))
temp2<-as.data.frame(as.matrix(sample_data(datITSr2)[,"Sample_name"]))
sort(as.numeric(temp2$Sample_name))
cbind(sort(as.numeric(as.character(temp$Sample_name))),sort(as.numeric(as.character(sort(temp2$Sample_name)))))
#temp3<-as.data.frame(as.matrix(sample_data(datBacr2)[,"Sample_name"]))

datEukr3<-subset_samples(datEukr2,Sample_name!="126")
datITSr3<-subset_samples(datITSr2,Sample_name!="81")

#calculated relative abundance of fungi
temp<-otu_table(datEukr3)
temp2<-aggregate.data.frame(temp,by=list(labels=labelsEuk),sum)
rownames(temp2)<-temp2$labels;temp2$labels<-NULL
fungirelabun<-t(temp2["Fungi",])

#remove fungi from euks
datEukr4<-subset_taxa(datEukr3,labels!="Fungi")#labels has no NAs so this works

datEukr4otu<-cbind(sample_data(datEukr4),t(otu_table(datEukr4)))
head(datEukr4otu)[,1:40]
datITSr3otu<-cbind(sample_data(datITSr3),t(otu_table(datITSr3)))
head(datITSr3otu)[,1:40]

#make sample numeric
datEukr4otu$Sample_name<-as.numeric(as.character(datEukr4otu$Sample_name))
datITSr3otu$Sample_name<-as.numeric(as.character(datITSr3otu$Sample_name))

#test if there is denovo name overlap, yes a ton, so relabel the denovos here and in the labels files
nameseuk<-names(datEukr4otu[,-c(1:31)])
namesits<-names(datITSr3otu[,-c(1:31)])
length(nameseuk)
length(namesits)
length(union(nameseuk,namesits))
intersect(nameseuk,namesits)

nameseuk2 <- sub("^", "e", nameseuk )
namesits2 <- sub("^", "i", namesits )

names(datEukr4otu)[-c(1:31)]<-nameseuk2
names(datITSr3otu)[-c(1:31)]<-namesits2

labelsEuk2<-labelsEuk
rownames(labelsEuk2)<-sub("^", "e",rownames(labelsEuk2))
labelsITS2<-labelsITS
rownames(labelsITS2)<-sub("^", "i",rownames(labelsITS2))



#sort things in numeric order by sample
fungirelabun2<-fungirelabun[order(datEukr4otu$Sample_name)]
names(fungirelabun2)<-rownames(fungirelabun)[order(datEukr4otu$Sample_name)]
datEukr5otu<-datEukr4otu[order(datEukr4otu$Sample_name),]
datITSr4otu<-datITSr3otu[order(datITSr3otu$Sample_name),]

#multiply ITS by fungi rel abun
datITSr4otusp<-datITSr4otu[,-c(1:31)]
datITSr4otusp2<-datITSr4otusp*fungirelabun2
rowSums(datITSr4otusp2)#this appears to have worked, the rowSumes (tot rel abun for each sample are the same as the fungirelabun2 file)

#add ITS to the right of the Euk data
datEukr6otu<-cbind(datEukr5otu,datITSr4otusp2)
dim(datEukr6otu)

#check
idenovo44940 in plot 5 is 0.00189401
temp<-otu_table(datITSr2)
head(temp)[,1:35]
temp[which(rownames(temp)=="denovo44940"),]
0.008008008, then scale *0.23651452 (from fungirelabun), checks


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









