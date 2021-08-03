# simulated examples in:
# "A One-dimension Genome Scan (1-DGS) Using the Mixed Model Approach Reveals QTL-by-genetic-background Interaction in Diallel and Nested Association Mapping designs"
# Wenhao Li & Martin Boer (2021)

library(asreml)
library(asremlPlus)
library(psych)
library(dplyr)
library(tidyr)
library(ggplot2)
source('R/functions.r')

dir=c('data/example_simDiallel')
# dir=c('data/example_simNAM')
progenyfile<-read.csv(paste0(dir,'/progenyfile.csv'))
progenyfile$pop<-factor(progenyfile$pop)
parGeno<-read.csv(paste0(dir,'/parGeno.csv'),stringsAsFactors = F)
progHap<-read.csv(paste0(dir,'/progHap.csv'),stringsAsFactors = F)
map<-read.csv(paste0(dir,'/map.csv'),stringsAsFactors = F)
progIBS<-read.csv(paste0(dir,'/progIBS.csv'),stringsAsFactors = F)

#________________________________________________________________________________________________________
# check some specific examples-----
# simulate epistasis e13=0.3 under a=0.3 or a=0 to check the QTL profile pattern and  
# validate the simulated QTL by model simQTL x Genome postion matrix
#_______________________________________________________________________________________________________________________________________________

##simulation
QTLnames=c('M57','M266','M440') #select three positions for analysis
parGeno %>% filter(marker %in% QTLnames)# parent genotypes at the three simQTL positions
a=0.3# additive effect for the 3 simQTLs
     # a=0# no additive effect for the 3 simQTLs
e13=rep(c(0.3),each=1)#add-add interaction effects between 1st and 3rd QTLs, one replication
parammat<- t(data.frame(a1=a,a2=a,a3=a,e12=0,e13=e13,e23=0,e123=0))

set.seed(1357)
allprogPheno<-simPheno(progIBS,QTLnames,parammat) #simulate phenotype
##analysis
threshold=-log10(0.05/nrow(map)) # Bonferroni-adjusted threshold
str.residual='hete' #'hete' --> family-specific（heterogeneous）residual structure 
                   #'homo' --> homogeneous redisudal structure
MQM=TRUE # TRUE --> multi-QTL model with cofactors (multiple rounds of genome scans)
         # FALSE --> single QTL model (one round of genome scan)
QTLwindow=20 # QTL window size for selecting cofactors (only works by setting MQM=TRUE )
trait.name='simY'

allMapRes<-list()
allprogPheno$set1$param #check paramters of this set
progPheno<-allprogPheno$set1$progPheno# the simY from this set
progPheno$pop<- factor(progPheno$pop)

## ref mod: IBD.MQM_F (Li et al., 2021)
allMapRes$refMod<-selQTL.IBD(parGeno,progHap,progPheno,map,trait.name,
                             QTLwindow,
                             threshold,
                             str.residual,
                             MQM)
## 1DGS mod
allMapRes$episMod<-selEpisQTL.IBD(parGeno,progHap,progPheno,map,trait.name,
                                  QTLwindow,
                                  threshold,
                                  str.residual,
                                  AIC.threshold=0)# by defalt 0 --> select the model with the lower AIC
                                                 # this value can be determined by permutation study with no-epistatic effect

plotEpisvsRef(allMapRes)

#validate simulated interacted QTL in biparent family
selpop<-'AxB' # can only be tested in the family where both simQTL1 and simQTL3 segregate
source('R/ValidateEpisQTL_BPP.R')
QxGplot

#________________________________________________________________________________________________________
# run 500 times for a grid of interacted effects-----
# only extract simQTL positions for testing in ref mod and 1DGS mod to save time
#_______________________________________________________________________________________________________________________________________________

QTLnames=c('M57','M266','M440') #select three positions for analysis
a=0.3# additive effect for the 3 simQTLs
    # a=0 # no additive effect for the 3 simQTLs

rep.n=500
e13=rep(c(0,0.1,0.2,0.3),each=rep.n)#add-add interaction effects between 1st and 3rd QTLs, 500 rep for each
parammat<- t(data.frame(a1=a,a2=a,a3=a,e12=0,e13=e13,e23=0,e123=0))

set.seed(1234)
allprogPheno<-simPheno(progIBS,QTLnames,parammat)
source('R/testEpisQTL_MPP.r')
sumBar


