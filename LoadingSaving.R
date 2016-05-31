#Loading/saving/packages needed

setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/NWT_MovingUphill2")
setwd("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata")

save.image("~/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill2_Workspace_Analysis.Rdata") 
load("/Users/farrer/Dropbox/EmilyComputerBackup/Documents/Niwot_King/Figures&Stats/kingdata/MovingUphill2_Workspace_Analysis.Rdata")


library(phyloseq)
#packageVersion("phyloseq")
library(foreach)
library(doParallel)


library(vegan)
library(reshape)
library(plotrix)
library(foreach)
library(doParallel)
library(Kendall)
library(tidyr)
library(grid)
library(data.table)
#library(BiodiversityR) #this requires X11 and takes a while to load, you need to close the window that it opens in rcommander
library(picante)
library(ggplot2)
library(NetIndices)
library(igraph)
library(fdrtool)
