
#sacn the two family-specific QTLs----

parents = names(parGeno)[-c(1:3)]
nloc = nrow(map)
megaM=design.HapM(progHap,parents,1:nloc)
allpop<-as.character(unique(progPheno$pop))


QTLcandidates=c(57,266,440)


for (c in 1:length(QTLcandidates)) {
  q=QTLcandidates[c]
  scanMat_q<-discon.designMat(sel_pos=(q),parGeno,progPheno,megaM)
  par1col<-seq(1,ncol(scanMat_q),by=2)
  par2col<-seq(1,ncol(scanMat_q),by=2)+1
  
  for (i in 1:length(par1col)) {
    z_q<-scanMat_q[,par1col[i]]-scanMat_q[,par2col[i]]
    z_sel<-scanMat_q[,par1col[i]]+scanMat_q[,par2col[i]]
    selRowM<-which(z_sel!=0)
    if (i ==1){
      newMat_q<-z_q[selRowM]
    } else {
      newMat_q<-c(newMat_q,z_q[selRowM])
    }
  }
  
  newMat_q<-as.matrix(newMat_q,ncol=1)
  newMat_q<-as.data.frame(newMat_q)
  names(newMat_q)<-paste0('simQTL',c)
  
  if (c==1) {
    allMat_q<-newMat_q
  } else {
    allMat_q<-cbind(allMat_q,newMat_q)
  }
  
}

for (c in 1:length(QTLcandidates)) {
  q=QTLcandidates[c]
  print(paste0('scan for simQTL', c,' x Genome interaction in family ', selpop))
  scanMat_q<-discon.designMat(sel_pos=(q),parGeno,progPheno,megaM)
  par1col<-seq(1,ncol(scanMat_q),by=2)
  par2col<-seq(1,ncol(scanMat_q),by=2)+1
  
  nloc<-ncol(progIBS)-1
  for (m in 1:nloc) {
    
    scanMat_m<-discon.designMat(sel_pos=m,parGeno,progPheno,megaM)
    for (i in 1:length(par1col)) {
      z_m<-scanMat_m[,par1col[i]]-scanMat_m[,par2col[i]]
      z_sel<-scanMat_m[,par1col[i]]+scanMat_m[,par2col[i]]
      selRow<-which(z_sel!=0)
      if (i ==1){
        newMat_m<-z_m[selRow]
      } else {
        newMat_m<-c(newMat_m,z_m[selRow])
      }
      
    }
    
    newMat_m<-allMat_q[,c]*newMat_m
    newMat_m<-as.matrix(newMat_m,ncol=1)
    if (m==1){
      allnewMat_m<-newMat_m
    } else {
      allnewMat_m<-cbind(allnewMat_m,newMat_m)
    }
    
  }
  
  allnewMat_m<-as.data.frame(allnewMat_m)
  names(allnewMat_m)<-paste0('QxM',1:nloc)
  
  allxmat<-cbind(allMat_q,allnewMat_m)
  data<-cbind(progPheno,allxmat)  
  data<-as.data.frame(data)
  
  # for (selpop in allpop){
  # print(selpop)
  popdata<-data%>% filter(pop==c(selpop))
  
  #model fit
  minlog10p<-c()
  for (m in 1:nloc){
    
    fix<-as.formula(paste0(trait.name, "~",paste0('simQTL1+simQTL2+simQTL3+QxM',m)))
    
    rcov<-as.formula(paste0("~idv(units)"))
    
    obj = asreml(fixed = fix, 
                 rcov = rcov, 
                 data = popdata, trace =FALSE)
    minlog10p[m]<--log10(anova(obj)[5,4])
    
  }
  res_QxGenome<-map %>%mutate(minlog10p_QxGenome=minlog10p,QTLcandidates=paste0('simQTL',c,' x Genome'),pop=selpop)
  
  if (q==QTLcandidates[1]) {
    allres_QxGenome<-res_QxGenome
  }else {
    allres_QxGenome<- rbind(allres_QxGenome,res_QxGenome)
  }
  # }
  
}

head(allres_QxGenome)

nchr = nlevels(factor(map$chr))
cmin = tapply(map$pos, map$chr, min)
cmax = tapply(map$pos, map$chr, max)
range = (cmax - cmin) + 8
sumrange = cumsum(range)
start_chr = c(0, sumrange[1:(nchr - 1)])
ver_line = sumrange[1:(nchr - 1)]
xtics = (start_chr + sumrange) * 0.5
x = start_chr[map$chr] + (map$pos - cmin[map$chr]) + 4
allres_QxGenome<-allres_QxGenome %>% mutate(x=rep(x,nrow(allres_QxGenome)/length(x)))

QxGplot<-ggplot(data=allres_QxGenome,aes(x=x,y=minlog10p_QxGenome,group=QTLcandidates))+
  geom_line(aes(linetype=QTLcandidates),size=0.7)+
  geom_point(data=allres_QxGenome %>% filter(name %in% QTLnames),aes(y=0,x=x),pch=4,size=3,color='red')+
  # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"),alpha=0.8)+
  # scale_fill_manual(values = alpha(c("#999999", "#E69F00", "#56B4E9"), 0.8))+
  geom_hline(yintercept=threshold,color='red',size=0.5)+
  geom_vline(xintercept=ver_line,color='grey10',lty=3,size=1)+
  # facet_grid(QTLcandidates~.)+
  scale_x_discrete(limits=tapply(x, map$chr, mean),
                   labels=1:max(allres_QxGenome$chr))+
  ggtitle(paste0('simQTL x Genome scan in family ',selpop))+
  xlab('Chromosome')+
  ylab('minlog10p')+
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.ticks.x = element_blank(),
    # legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.border = element_rect(fill = NA, color = "black",size = 0.5, linetype = "solid")
  )

