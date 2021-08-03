#example to transfer the format required for tha analysis in:
#"A One-dimension Genome Scan Using the Mixed Model Approach Reveals QTL-by-genetic-background Interaction in Diallel and Nested Association Mapping designs"
# Wenhao Li & Martin Boer (2021)

##load orinial files
library(asreml)
library(asremlPlus)
library(psych)
library(dplyr)
library(tidyr)
library(ggplot2)
source('functions.r')

dir='data/example_tomatoDiallel/originalfiles'
pedInfofile='tomatoDiallel_pedinfor.csv'
magicfile='tomatoDiallel_magicReconstruct_HaploProb.csv'
phenofile='phenofile.csv'
sel.dist=0.1

# # read data

pedinfo<-read.table(paste0(dir,'/',pedInfofile),
                    header = T,
                    sep = ",",
                    skip = 1,
                    stringsAsFactors = F)
generation0<-pedinfo %>% filter(pedinfo$Generation==0)
npar<-nrow(generation0)

magicHap<-read.csv(paste0(dir,'/',magicfile),
                   header = F,
                   stringsAsFactors = F)
map<-data.frame(name=as.character(magicHap[2,-1]),
                chr=as.numeric(magicHap[3,-1]),
                pos=as.numeric(magicHap[4,-1]))
parGeno<-as.data.frame(t(magicHap[2:(which(magicHap[,2]=='logl')-1),-1]))
names(parGeno)<-magicHap[2:(which(magicHap[,2]=='logl')-1),1]
progHap<-magicHap[-c(1:(which(magicHap[,2]=='haploprob')+3)),]
colnames(progHap)<-c('prog_hap',as.character(map$name))

progPheno<-read.csv(paste0(dir,'/',phenofile))
progPheno$pop<-factor(progPheno$pop)
## select markers for mapping
sel.map<-select.markers(map,sel.dist)
sel.progHap<- progHap[,c(1,which(names(progHap) %in% sel.map$name))]
sel.parGeno<-parGeno[which(parGeno$marker %in% sel.map$name),]
#
write.csv(sel.parGeno,'data/example_tomatoDiallel/parGeno.csv',row.names = F)
write.csv(sel.progHap,'data/example_tomatoDiallel/progHapsel.csv',row.names = F)
write.csv(sel.map,'data/example_tomatoDiallel/map.csv',row.names = F)