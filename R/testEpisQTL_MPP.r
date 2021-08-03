parents = names(parGeno)[-c(1:3)]
nloc = nrow(map)
megaM=design.HapM(progHap,parents,1:nloc)
threshold=-log10(0.01)
sel.markers<-c(57,266,440)
AIC_threshold<-0
for ( j in 1:length(names(allprogPheno))) {
  print(paste0('for set:',j))
  allprogPheno[[paste0('set',j)]]$param #check paramters of this set
  progPheno<-allprogPheno[[paste0('set',j)]]$progPheno# the simY from this set
  progPheno$pop<- factor(progPheno$pop)
  
  #refMod————————————————————————————————————————————————————————————————————————————————————————----
  selcof=NULL
  minlog10p<-c()
  refRes<-data.frame(simQTL=paste0('simQTL',1:3),sel.markers=sel.markers,
                     cofactor=FALSE)
  
  for (r in 1:4) {
    print(r)
    for (m in 1:length(sel.markers)) {
      df2<-refRes[-m,]
      selcof<-df2$sel.markers[which(df2$cofactor)]
      if (length(selcof)==0) {
        selcof<-NULL
      }
      obj1<-multiHAPmodel(parGeno,megaM,progPheno, trait.name='simY', scan_pos=sel.markers[m],cof=selcof,str.residual='hete',null.model=F)
      obj0<-multiHAPmodel(parGeno,megaM,progPheno, trait.name='simY', scan_pos=sel.markers[m],cof=selcof,str.residual='hete',null.model=T)
      
      dev0 = -2.0*obj0$asreml.obj$loglik
      dev1 = -2.0*obj1$asreml.obj$loglik
      dev = dev0 - dev1
      minlog10p[m] = -log10(0.5*pchisq(dev,1,lower.tail=FALSE))
      # print(paste0('select cof:', selcof, 'for',sel.markers[m]))
    }
    # print(paste0('minlog10p: ',minlog10p,'in rep ',r))
    refRes$minlog10p<-minlog10p
    if (length(which(refRes$cofactor))==0) {
      if (max(minlog10p)<threshold) break
      cof.ind<-which(minlog10p==max(minlog10p))
      refRes$cofactor[cof.ind]<-TRUE
    } else {
      if (length(which(refRes$cofactor))==3) break
      if (max(minlog10p[-which(refRes$cofactor)])<threshold) break
      cof<-which(minlog10p==max(minlog10p[-which(refRes$cofactor)]))
      refRes$cofactor[cof]<-TRUE
    }
  }
  refRes$set=j
  if (j==1) {
    allrefRes<-refRes
  } else {
    allrefRes<-rbind(allrefRes,refRes)
  }
  
  #1DGS Mod————————————————————————————————————————————————————————————————————————————————————————----
  selcof=NULL
  minlog10p<-c()
  type<-c()
  AIC.diff<-c()
  episRes<-data.frame(simQTL=paste0('simQTL',1:3),sel.markers=sel.markers,
                      cofactor=FALSE,cofactor.defin='unknown')
  for (r in 1:4) {
    print(r)
    for (m in 1:length(sel.markers)) {
      df2<-episRes[-m,]
      selcof<-df2$sel.markers[which(df2$cofactor)]
      definselCof<-df2$cofactor.defin[which(df2$cofactor)]
      if (length(selcof)==0) {
        selcof<-NULL
        definselCof<-NULL
      }
      
      obj_parent<-generalEpis.model(trait.name='simY',scan_pos=sel.markers[m],definScan='parent',
                                    cof=selcof,definCof=definselCof,
                                    parGeno,
                                    progPheno,
                                    megaM,
                                    str.residual='hete',
                                    null.model=FALSE)
      obj_family<-generalEpis.model(trait.name='simY',scan_pos=sel.markers[m],definScan='family',
                                    cof=selcof,definCof=definselCof,
                                    parGeno,
                                    progPheno,
                                    megaM,
                                    str.residual='hete',
                                    null.model=FALSE)
      
      obj_null<-generalEpis.model(trait.name='simY',scan_pos=sel.markers[m],definScan='family',
                                  cof=selcof,definCof=definselCof,
                                  parGeno,
                                  progPheno,
                                  megaM,
                                  str.residual='hete',
                                  null.model=TRUE)
      
      dev_main_par<--2*(obj_null$asreml.obj$loglik-obj_parent$asreml.obj$loglik)
      (minlog10p_main_par<--log10(0.5*pchisq(dev_main_par,1,lower.tail=FALSE)))
      
      dev_main_fam<--2*(obj_null$asreml.obj$loglik-obj_family$asreml.obj$loglik)
      (minlog10p_main_fam<--log10(0.5*pchisq(dev_main_fam,1,lower.tail=FALSE)))
      
      
      IC.family<- infoCriteria.asreml(obj_family$asreml.obj)
      IC.parent<- infoCriteria.asreml(obj_parent$asreml.obj)
      (AIC_epis<-(IC.parent$AIC-IC.family$AIC))
      AIC.diff[m]<-AIC_epis
      if(AIC_epis< AIC_threshold) {
        minlog10p[m]<-minlog10p_main_par
        type[m]<-'parent'
      } else {
        minlog10p[m]<-minlog10p_main_fam
        type[m]<-'family'
      }
      
    }
    episRes$AIC.diff<-AIC.diff
    episRes$minlog10p<-minlog10p
    episRes$cofactor.defin<-type
    if (length(which(episRes$cofactor))==0) {
      if (max(minlog10p)<threshold) break
      cof.ind<-which(minlog10p==max(minlog10p))
      episRes$cofactor[cof.ind]<-TRUE
    } else {
      if (length(which(episRes$cofactor))==3) break
      if (max(minlog10p[-which(episRes$cofactor)])<threshold) break
      cof<-which(minlog10p==max(minlog10p[-which(episRes$cofactor)]))
      episRes$cofactor[cof]<-TRUE
    }
  }
  
  episRes$set=j
  if (j==1) {
    allepisRes<-episRes
  } else {
    allepisRes<-rbind(allepisRes,episRes)
  }
}


#summary----

head(allrefRes)
head(allepisRes)
threshold
maxepis=0
rep.n=500
set.n<-length(names(allprogPheno))/rep.n

for (s in 1:set.n) {
  
  set1refRes<-allrefRes %>% filter(set %in% ((s-1)*rep.n+1):(s*rep.n)) 
  set1Res<-allepisRes %>% filter(set %in% ((s-1)*rep.n+1):(s*rep.n)) %>% 
    mutate(ref.minlog10p=set1refRes$minlog10p,epis=paste0('episSet',s))
  
  for (q in 1:3){
    #by episMOd
    episModQTL<-set1Res%>% filter(simQTL==paste0('simQTL',q),minlog10p>threshold) 
    
    episModQTL_f<-episModQTL%>% filter(cofactor.defin=='family')
    nrow(episModQTL_f)
    episModQTL_p<-episModQTL%>% filter(cofactor.defin=='parent')
    nrow(episModQTL_p)
    #by refMod
    refModQTL<-set1Res%>% filter(simQTL==paste0('simQTL',q),ref.minlog10p>threshold) 
    nrow(refModQTL)
    
    nQTL<-data.frame(simQTL=paste0('simQTL',q),mod=c('1DGS','1DGS','IBD.MQM_F'),nQTL=c(nrow(episModQTL_f),
                                                                                       nrow(episModQTL_p),
                                                                                       nrow(refModQTL)),
                     defin=c('family','parent','ref'),epis=paste0('episSet',s))
    
    if (q==1){
      nQTLsum<-nQTL
    } else {
      nQTLsum<-rbind(nQTLsum,nQTL)
    }
  }
  
  if (s==1){
    Res<-set1Res
    allnQTLsum<-nQTLsum
  } else {
    Res<-rbind(Res,set1Res)
    allnQTLsum<-rbind(allnQTLsum,nQTLsum)
  }
  
}

ggplot(data =Res,aes(x=ref.minlog10p,y=minlog10p,color=cofactor.defin) )+
  geom_point(size=0.8,alpha=0.7)+
  geom_vline(xintercept=threshold,color='grey10',lty=2,size=0.5)+
  geom_hline(yintercept=threshold,color='grey10',lty=2,size=0.5)+
  facet_grid(epis~simQTL)+
  ggtitle(paste0('maxEpis:',maxepis,'+QTLthreshold (-log10p):',threshold,' + AICthreshold:', AIC_threshold))+
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    plot.background = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.position = "top",
    strip.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black",size = 0.5, linetype = "solid"))


allnQTLsum$nQTL<-allnQTLsum$nQTL/rep.n #transfer to percentages
sumBar<-ggplot(allnQTLsum, aes(x=epis,y=nQTL,fill=defin)) +
  geom_bar(stat = "identity",width = 0.3)+
  scale_fill_manual(values=c("coral2", "cyan3", "grey60"))+
  ylim(0,1)+
  xlab('detection power')+
  facet_grid(mod~simQTL)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "top",
    strip.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "white",size = 0.5, linetype = "solid"))




