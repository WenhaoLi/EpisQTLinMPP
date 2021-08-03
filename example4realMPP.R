# Real examples for analysis in:
# "A One-dimension Genome Scan Using the Mixed Model Approach Reveals QTL-by-genetic-background Interaction in Diallel and Nested Association Mapping designs"
# Wenhao Li & Martin Boer (2021)

library(asreml)
library(asremlPlus)
library(psych)
library(dplyr)
library(tidyr)
library(ggplot2)
source('R/functions.r')

# read data----
# dir='data/example_tomatoDiallel'
dir='data/example_maizeDiallel'# markers were selected at 5cM to reduce analysis time as an example
progPheno<-read.csv(paste0(dir,'/phenofile.csv'))
progPheno$pop<-factor(progPheno$pop)
parGeno<-read.csv(paste0(dir,'/parGeno.csv'),stringsAsFactors = F)
progHap<-read.csv(paste0(dir,'/progHap.csv'),stringsAsFactors = F)
map<-read.csv(paste0(dir,'/map.csv'),stringsAsFactors = F)

#analysis----
threshold=-log10(0.05/nrow(map)) # Bonferroni-adjusted threshold
str.residual='hete' #'hete' --> family-specific（heterogeneous）residual structure 
                    #'homo' --> homogeneous redisudal structure
MQM=TRUE # TRUE --> multi-QTL model with cofactors (multiple rounds of genome scans)
         # FALSE --> single QTL model (one round of genome scan)

QTLwindow=20 # QTL window size for selecting cofactors (only works by setting MQM=TRUE )

names(progPheno) 
trait.name='PHS' #obtain the trait name for analysis

allMapRes<-list()
## Reference model
allMapRes[[trait.name]]$refMod<-selQTL.IBD(parGeno,
                                           progHap,
                                           progPheno,
                                           map,
                                           trait.name,
                                           QTLwindow,
                                           threshold,
                                           str.residual,
                                           MQM)
## map epistitic QTL
allMapRes[[trait.name]]$episMod<-selEpisQTL.IBD(parGeno,
                                                progHap,
                                                progPheno,
                                                map,
                                                trait.name,
                                                QTLwindow,
                                                threshold,
                                                str.residual,
                                                MQM,
                                                AIC.threshold=0)# by defalt 0 --> select the model with the lower AIC
## Plot
plotEpisvsRef(allMapRes[[trait.name]])

#calculate the R squre explained by random effect (Piepho 2019)----
## R^2 by random QTLs by refMod model
candQTL_refMod<-list(cofactors=allMapRes[[trait.name]]$refMod[[1]],
                     definCof=rep('parent',length(allMapRes[[trait.name]]$refMod[[1]])))
candQTL_refMod
round(calcR_sq_random(candQTL_refMod,trait.name,parGeno,progPheno,progHap,map),2)

## R^2 by random QTLs by episMod model
candQTL_episMod<-allMapRes[[trait.name]]$episMod[[1]]
candQTL_episMod
round(calcR_sq_random(candQTL_episMod,trait.name,parGeno,progPheno,progHap,map),2)
