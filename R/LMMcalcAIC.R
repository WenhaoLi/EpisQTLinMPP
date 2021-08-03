# AIC calculation in linear mixed model from asremlR(v3) object
# AIC of ASREML-r(v4) object can be directly extracted
# Wenhao 2021.7.19


dir='data/example_maizeDiallel'
progPheno<-read.csv(paste0(dir,'/phenofile.csv'))
progPheno$pop<-factor(progPheno$pop)
parGeno<-read.csv(paste0(dir,'/parGeno.csv'),stringsAsFactors = F)
progHap<-read.csv(paste0(dir,'/progHap.csv'),stringsAsFactors = F)
map<-read.csv(paste0(dir,'/map.csv'),stringsAsFactors = F)


parents = names(parGeno)[-c(1:3)]
nloc = nrow(map)
megaM=design.HapM(progHap,parents,1:nloc)

scan_pos<- 1 # fit the 1st marker as example
definScan<-c('parent','family')[1]

cof=NULL # no cofactor
definCof=NULL

pop = "pop"
geno = "geno"
sel_pos=c(scan_pos,cof)
npos = length(sel_pos)
npar=length(parents)

## design matrix
if (definScan=='family') {
  scanMat<-discon.designMat(scan_pos,parGeno,progPheno,megaM)
} 

if (definScan=='parent') {
  scanMat<-con.designMat(scan_pos,parGeno,megaM)
}

#data to fit wit asreml
data = cbind(progPheno, scanMat)

## R structure
rcov = as.formula(paste0("~at(pop):units"))

## fixed effect
fix = as.formula(paste0(trait.name, "~", pop))

## group for random effect
sel_defin<-c(definScan,definCof)

allMnames<-c()
for (q in 1:length(sel_pos)) {
  
  if (sel_defin[q]=='parent'){
    Mnames=paste0('M',sel_pos[q],'p_')
    
  } else {
    Mnames<-paste0('M',sel_pos[q],paste0('f',1:length(unique(progPheno$pop))),'_')
    
  }
  allMnames<-c(allMnames,Mnames)
}

formQTLs<- paste0("grp(", allMnames, ")", collapse = "+")
group=list()
for (n in 1:length(allMnames)){
  group[[allMnames[n]]]<-grep(paste0(allMnames[n]),names(data))
}

##random term
random=as.formula(paste0('~ ',formQTLs))

##fit mixed model
obj = asreml(fixed = fix, 
               rcov = rcov, 
               random=random,
               group=group[],
               data = data, trace =FALSE)
assign(paste0('obj_',definScan),obj)
summary(obj)

## use asremlPlus (Chris Brien) 
(IC<-infoCriteria.asreml(obj,IClikelihood = "REML"))
IC$AIC

## (Boer et al. 2020)
2*(IC$varDF)-2*IC$loglik

IC$AIC==2*(IC$varDF)-2*IC$loglik # they are equal






