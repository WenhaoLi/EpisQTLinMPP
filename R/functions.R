# Functions required for:
# "A One-dimension Genome Scan Using the Mixed Model Approach Reveals QTL-by-genetic-background Interaction in Diallel and Nested Association Mapping designs" 
# Wenhao Li & Martin Boer (2021)

# Required functions-----

#the function transfers the IBDs files into the formated of design matrix
design.HapM<-function(progHap,parents,sel_pos) {
  Hap_mat<-progHap[,-1]
  for(i in 1:length(sel_pos)) {
    HapPro<-as.numeric(Hap_mat[,sel_pos[i]])
    M1<-matrix(HapPro,ncol = length(parents),byrow=T)
    M1<-as.data.frame(M1)
    colnames(M1)<-paste0('haplotype',1:length(parents),'M',sel_pos[i])
    rownames(M1)<- unique(substr(progHap$prog_hap,1,regexpr('_haplotype',progHap$prog_hap)-1))
  
    if(i==1) {
      M=M1
    } else {
      M<-cbind(M,M1)
    }
  }
  return(M)
}

#the function plots QTL profile during genome scans
plotQTLscan = function(data, threshold, cofactors, trait.name)
{
  title = paste('QTL-profile for trait ', trait.name)
  if (length(cofactors)==0) {
    title = paste0(title,', no cofactors')
  } else if (length(cofactors) == 1) {
    title = paste0(title,', one cofactor')
  } else {
    title = paste0(title,', ',length(cofactors),' cofactors')
  }
  
  nchr = nlevels(factor(data$chr))
  cmin = tapply(data$pos, data$chr, min)
  cmax = tapply(data$pos, data$chr, max)
  range = (cmax-cmin)+8.0
  sumrange = cumsum(range)
  start_chr = c(0,sumrange[1:(nchr-1)])
  ver_line = sumrange[1:(nchr-1)]
  xtics = (start_chr+sumrange)*0.5
  x = start_chr[data$chr] + (data$pos-cmin[data$chr])+4
  
  chromF = as.factor(data$chr)
  col_chr = c("black","red","green","blue","cyan1","purple","gold","brown","chartreuse","black",
              "black","black")
  col=col_chr[chromF]
  max_y = 5.0*ceiling(max(data$minlog10p)/5.0)
  plot(x=x,y=data$minlog10p,ylim=c(0,max_y),col='blue',yaxs="i",cex = 0.55,type='l',
       axes=FALSE,xpd=FALSE,pch=16,ylab='-log10(P)',xlab='Chromosomes',main=title)
  # defines the end points of the chromosomes:
  abline(v=ver_line,col='red', lwd = 2.0,lty=3)
  abline(h=threshold,col='red',lwd = 0.7)
  axis(1,at=xtics,labels=c(1:nchr))
  axis(2)
  box()
}

#the function visually compares the refMod and episMod
plotEpisvsRef<-function(MapRes){
  data=MapRes$refMod[[2]]
  nchr = nlevels(factor(data$chr))
  cmin = tapply(data$pos, data$chr, min)
  cmax = tapply(data$pos, data$chr, max)
  range = (cmax - cmin) + 8
  sumrange = cumsum(range)
  start_chr = c(0, sumrange[1:(nchr - 1)])
  ver_line = sumrange[1:(nchr - 1)]
  xtics = (start_chr + sumrange) * 0.5
  x = start_chr[data$chr] + (data$pos - cmin[data$chr]) + 4
  
  RefRes<- MapRes$refMod[[2]]
  EpisRes<- MapRes$episMod[[2]]
  cofactors<- MapRes$episMod[[1]]
  
  p<-ggplot(data=RefRes, aes(x=x,y=minlog10p))+
    geom_area(fill="grey60",alpha=1)+
    geom_point(data=EpisRes,aes(x=x,y=minlog10p,color=QTLdefin))+
    # scale_color_manual(values=c("lightcoral", "skyblue"))
    geom_hline(yintercept=threshold,color='red',size=0.5)+
    geom_vline(xintercept=ver_line,color='grey10',lty=3,size=1)+
    ggtitle(paste('Epistatic QTL mapping  \n vs. \n reference model (grey background)'))+
    scale_x_discrete(limits=tapply(x, RefRes$chr, mean),
                     labels=as.character(unique(RefRes$chr)))+
    xlab('Chromosome')+
    ylab('minlog10p')+
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      strip.background = element_blank(),
      panel.border = element_rect(fill = NA, color = "black",size = 0.5, linetype = "solid"))
	  
  print(p)
}

#the function selects marker at a certain distance from a genetic map
select.markers<- function(test.map,min.dist) {
  final.map=NULL
  if (min.dist==0) {
    final.map=test.map
  }
  if(!min.dist==0) {
    all.chr<-unique(test.map$chr)
    for( c in all.chr) {
      chr.map=test.map[which(test.map$chr==c),]
      m=1
      chr.pos<-chr.map$pos
      sel.markers<-NULL
      for (i in 1:length(chr.pos)) {
        sel.markers[i]=m
        next.m<-which(!(chr.pos-chr.pos[m])<min.dist)[1]
        if(is.na(next.m)) break
        m=next.m
      }
      sel.chr.map<-chr.map[sel.markers,]
  
      if(c==all.chr[1]) {
        final.map<-sel.chr.map
      } else {
        final.map<-rbind(final.map,sel.chr.map)
      }
    }
  }
  
  return(final.map)
}

# the fucntions select cofactors outside of a certain QTLwindow size
calc.distance = function(markers,m1,m2)
{
  dist = ifelse(markers[m1,]$chr==markers[m2,]$chr,
                abs(markers[m1,]$pos - markers[m2,]$pos), 1.0e5)
  dist
}
select.cofactors = function(markers, m, cofactors, QTLwindow)
{
  if (length(cofactors)==0) return(NULL)
  min.dist = 0.5*QTLwindow
  dist = sapply(cofactors,function(x) calc.distance(markers,m,x))
  ord = order(dist)
  if (dist[ord[1]] > min.dist)
  {
    return(cofactors)
  }
  cofactors[-c(ord[1])]
}

#the function geneartes the indications in the design matrix for randon effects in one group
generate.group = function(ncol_pheno,nparents,names)
{
  group = list()
  for (c in 1:length(names))
  {
    group[[names[c]]] = seq(from=ncol_pheno+1+nparents*(c-1),length=nparents) 
  }
  group  
}
#the function simulate phenotypes
simPheno<-function(progIBS,QTLnames,parammat){
  allprogPheno<-list()
  for (i in 1:ncol(parammat)) {
    (a1=parammat[1,i])
    (a2=parammat[2,i])
    (a3=parammat[3,i])
    (e12=parammat[4,i])
    (e13=parammat[5,i])
    (e23=parammat[6,i])
    (e123=parammat[7,i])
    
    qtlIBS<-progIBS %>% select(QTLnames)
    apply(qtlIBS, 2, table)
    # construct the beta matrix of genetic indicators
    j1<-qtlIBS[,1]
    j2<-qtlIBS[,2]
    j3<-qtlIBS[,3]
    
    x<-data.frame(z1=j1-1,
                  z2=j2-1,
                  z3=j3-1,
                  I12=(j1-1)*(j2-1),
                  I13=(j1-1)*(j3-1),
                  I23=(j2-1)*(j3-1),
                  I123=(j1-1)*(j2-1)*(j3-1))
    
    beta=matrix(c(a1,a2,a3,e12,e13,e23,e123),ncol=1)
    
    # pheno=x*beta+error
    
    simY=as.matrix(x) %*% beta + rnorm(nrow(x),0,0.6)
    
    progPheno<-data.frame(geno=progenyfile$geno,simY=simY,pop=progenyfile$pop)
    allprogPheno[[paste0('set',i)]]$param<-c(parammat[,i])
    allprogPheno[[paste0('set',i)]]$progPheno<-progPheno
  }
  return(allprogPheno)
}
# Functions for IBD-based multi-QTL (IBD.MQM) approach (as reference model in our study) ----

multiHAPmodel<-function (parGeno,megaM,progPheno, trait.name, scan_pos,cof,str.residual,null.model=F)
{
  pop = "pop"
  geno = "geno"
  
  sel_pos = (c(scan_pos,cof))
  npos=length(sel_pos)
  fix = as.formula(paste0(trait.name, "~", pop))
  
  if (str.residual=='homo') {
    rcov = as.formula("~idv(units)")
  }
  
  if (str.residual=='hete') {
    rcov = as.formula(paste0("~at(pop):units"))
  }

  Mnames = c(paste0("M", sel_pos))
  formQTLs = paste0("grp(", Mnames, ")", collapse = "+")
  
  parents = names(parGeno)[-c(1:3)]
  
  sel_col=c()
  for (s in sel_pos){
    selM_col<-seq((s-1)*length(parents)+1,length.out = length(parents))
    sel_col<-c(sel_col,selM_col)
  }

  M<-as.data.frame(megaM[,sel_col])
  M<- M %>% mutate(geno=rownames(M))
  data = merge(progPheno, M,by.x = geno,sort = FALSE)

  ncol_pheno = ncol(progPheno)
  npar = length(parents)
  data = data.frame(data)
  group = generate.group(ncol_pheno, npar, Mnames)
  
  if(null.model==FALSE) {
    random=as.formula(paste0('~ ',formQTLs))
    capture.output({
      oldw <- getOption("warn")
      options(warn = -1)
      obj = asreml(fixed = fix, 
                   random = random, 
                   rcov = rcov, 
                   data = data, 
                   group = group[], 
                   trace = FALSE)
      options(warn = oldw)
    })
  }
  if(null.model==TRUE) {
    if(length(cof)==0) {
      capture.output({
        oldw <- getOption("warn")
        options(warn = -1)
        obj = asreml(fixed = fix, 
                     rcov = rcov, 
                     data = data, 
                     trace = FALSE)
        options(warn = oldw)
      })
    } else {
      Cofnames = c(paste0("M", cof))
      forcofLs = paste0("grp(", Cofnames, ")", collapse = "+")
      random=as.formula(paste0('~ ',forcofLs))
      capture.output({
        oldw <- getOption("warn")
        options(warn = -1)
        obj = asreml(fixed = fix, 
                     random = random, 
                     rcov = rcov, 
                     data = data, 
                     group = group[], 
                     trace = FALSE)
        options(warn = oldw)
      })
    }
    
  }
  obj
  return(list(data=data,asreml.obj=obj))
}

CIMscan.HAPrandom<-function(parGeno,megaM, progPheno, trait.name,map,cof=NULL,QTLwindow=30,str.residual) {
  parents = names(parGeno)[-c(1:3)]
  nloc = nrow(map)
  sel_markers=c(1:nloc)
  Loci=map
  
  minlog10p = vector(length=nloc)
  QTL_region = rep(FALSE,nloc)
  effects = matrix(nrow=nloc,ncol=length(parents))
  colnames(effects) = parents
  
  for (i in 1:nloc)
  {
    m = sel_markers[i]
    sel_cof = select.cofactors(Loci, m, cof, QTLwindow)
    if (length(sel_cof)!=length(cof)) {
      QTL_region[i] = TRUE
    }
    
    obj1=multiHAPmodel(parGeno,megaM,progPheno, trait.name, scan_pos=m,cof=sel_cof,str.residual,null.model=F)
    obj0=multiHAPmodel(parGeno,megaM,progPheno, trait.name, scan_pos=m,cof=sel_cof,str.residual,null.model=T)
    
    obj1<-obj1$asreml.obj
    obj0<-obj0$asreml.obj

    dev0 = -2.0*obj0$loglik
    dev1 = -2.0*obj1$loglik
    dev = dev0 - dev1
    minlog10p[i] = -log10(0.5*pchisq(dev,1,lower.tail=FALSE))
    scan_posname<-paste0('grp(M',i,')')

    if (i %% 25 == 0) {
      print(i)
    }
  }
  
  df = data.frame(ndx=sel_markers, Loci[sel_markers,],minlog10p=minlog10p,eff=effects, QTLregion = QTL_region, trait=trait.name)
  return(df)
}

selQTL.IBD<-function(parGeno,
                     progHap,
                     progPheno,
                     map,
                     trait.name,
                     QTLwindow,
                     threshold,
                     str.residual,
                     MQM){
  
  parents = names(parGeno)[-c(1:3)]
  nloc = nrow(map)
  cofactors = NULL

  megaM=design.HapM(progHap,parents,1:nloc)
  for (i in 1:100) {
    print(paste0("QTL scan for trait ", trait.name,", ", length(cofactors), " cofactors"))
    df.scan = CIMscan.HAPrandom(parGeno,megaM, progPheno, trait.name,map,cof=cofactors,QTLwindow,str.residual)
    # (parGeno,progHap, progPheno, trait.name,map,cof=NULL,QTLwindow=30,str.residual)# df.scan = MP.CIMscan(MP.map,MP.IBD,MP.pheno,trait.name,sel_markers=c(1:nloc),cof=cofactors,QTL.random = TRUE,QTLwindow ,str.residual)
    plotQTLscan(df.scan,threshold=threshold,cofactors,trait.name=trait.name)
    
    df.scan.sel = filter(df.scan, QTLregion==FALSE)
    # remove NAs from minlog10p values
    df.scan.sel = filter(df.scan.sel, is.na(minlog10p) == FALSE)
    ord = order(df.scan.sel$minlog10p, decreasing = TRUE)
    max_value = df.scan.sel[ord[1],]$minlog10p
    
    if (max_value < threshold) {
      if (!is.null(cofactors)) cofactors <- sort(cofactors)
    }
    
    if (max_value < threshold) break
    
    if (MQM==F) break
    
    cofactors = c(cofactors, df.scan.sel[ord[1],]$ndx)
    if (length(cofactors)==10) break
    
  }
  
  results <- list(cofactors, df.scan)
  return(results)
}

#Functions for IBD-based one-dimensional genome scan approach for detecting epistatic QTLs ----
con.designMat<- function (sel_pos,parGeno,megaM) {
  
  parents = names(parGeno)[-c(1:3)]
  npar=length(parents)
  for (s in sel_pos){
    selM_col<-seq((s-1)*length(parents)+1,length.out = length(parents))
    M_scan<-as.data.frame(megaM[,selM_col])
    
    if (s==sel_pos[1]) {
      allM_scan<-M_scan
    } else {
      allM_scan<- cbind(allM_scan,M_scan)
    }
  }
  
  colnames(allM_scan)<-paste0(colnames(allM_scan),'p_')
  return(allM_scan)
}

discon.designMat<-function(sel_pos,parGeno,progPheno,megaM){
  parents = names(parGeno)[-c(1:3)]
  npar=length(parents)
  for (s in sel_pos){
    selM_col<-seq((s-1)*length(parents)+1,length.out = length(parents))
    M<-as.data.frame(megaM[,selM_col])
    fnames<- unique(progPheno$pop)
    newM<-matrix(0,ncol=2*length(fnames),nrow=nrow(progPheno))
    newM<- as.data.frame(newM)
    for (f in 1:length(fnames)) {
      frow<-which(progPheno$pop==fnames[f])
      frowM<-M[frow,]
      SelectM<-frowM[, colSums(frowM != 0) > 0]
      newM[frow,2*f-1]<-SelectM[,1]
      newM[frow,2*f]<-SelectM[,2]
      names(newM)[c(2*f-1,2*f)]<-paste0(names(SelectM),'f',f,'_')
    }
    
    head(newM)
    
    
    if (s==sel_pos[1]) {
      allnewM<-newM
    } else {
      allnewM<- cbind(allnewM,newM)
    }
    
  }
  
  return(allnewM)
}

generalEpis.model<- function(trait.name,
                             scan_pos=1,definScan=c('parent','family')[2],
                             cof=NULL,definCof=NULL,
                             parGeno,
                             progPheno,
                             megaM,
                             str.residual,
                             null.model=FALSE) {
  
  
  pop = "pop"
  geno = "geno"
  sel_pos=c(scan_pos,cof)
  npos = length(sel_pos)
  parents = names(parGeno)[-c(1:3)]
  npar=length(parents)
  
  ## design matrix
  if (definScan=='family') {
    scanMat<-discon.designMat(scan_pos,parGeno,progPheno,megaM)
  } else {
    scanMat<-con.designMat(scan_pos,parGeno,megaM)
  }
  
  if (length(cof)!=0){
    
    for (c in 1:length(cof)) {
      
      if (definCof[c]=='family') {
        
        cofMat<-discon.designMat(cof[c],parGeno,progPheno,megaM)
      } else {
        cofMat<-con.designMat(cof[c],parGeno,megaM)
      }
      
      if (c==1) {
        allcofMat<-cofMat
      } else {
        allcofMat<-cbind(allcofMat,cofMat)
      }
      
    }
    
    allxmat<-cbind(scanMat,allcofMat)
    
  } else {
    allxmat<-scanMat
  }
  
  # allxmat$geno<-rownames(allxmat)
  data = cbind(progPheno, allxmat)
  
  ## R structure
  if (str.residual=='homo') {
    rcov = as.formula("~idv(units)")
  }
  
  if (str.residual=='hete') {
    rcov = as.formula(paste0("~at(pop):units"))
  }
  
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
  
  if (null.model==FALSE) {
    random=as.formula(paste0('~ ',formQTLs))
    obj = asreml(fixed = fix, 
                 rcov = rcov, 
                 random=random,
                 group=group[],
                 data = data, trace =FALSE)
  }
  
  if (null.model==TRUE) {
    if(length(cof)==0) {
      obj = asreml(fixed = fix, 
                   rcov = rcov, 
                   data = data, trace =FALSE)
    } else {
      allMnames<-c()
      for (q in 1:length(cof)) {
        
        if (definCof[q]=='parent'){
          Mnames=c(paste0('M',cof[q],'p_'))
          allMnames<-c(allMnames,Mnames)
          
        } else {
          Mnames<-paste0('M',cof[q],paste0('f',1:length(unique(progPheno$pop))),'_')
          allMnames<-c(allMnames,Mnames)
        }
      }
      
      formQTLs<- paste0("grp(", allMnames, ")", collapse = "+")
      random=as.formula(paste0('~ ',formQTLs))
      obj = asreml(fixed = fix, 
                   rcov = rcov, 
                   random=random,
                   group=group[],
                   data = data, trace =FALSE)
    }
    
  }
  
  list.obj=list()
  list.obj[['data']]<-data
  list.obj[['asreml.obj']]<-obj
  return(list.obj)
  
}

CIMscan.genEpis<-function(trait.name='simY',
                          megaM,
                          parGeno,
                          progPheno,
                          map,
                          cof=NULL,
                          definCof,
                          QTLwindow=20,
                          str.residual='hete',
                          AIC.threshold=0){
  
  nloc = nrow(map)
  minlog10p_main_par=c()
  minlog10p_main_fam=c()
  AIC_epis=c()
  QTL_region = rep(FALSE,nloc)
  print(paste0('scan genome ', ' with ',length(cof),' cofactor(s)'))
  for(m in 1:nloc) {
    sel_cof = select.cofactors(map, m, cof, QTLwindow)
    sel_definCof=definCof[which(cof %in% sel_cof)]
    
    if (length(sel_cof)!=length(cof)) {
      QTL_region[m] = TRUE
    }
    
    obj_parent<-generalEpis.model(trait.name,
                                  scan_pos=m,definScan='parent',
                                  cof=sel_cof,definCof=sel_definCof,
                                  parGeno,
                                  progPheno,
                                  megaM,
                                  str.residual,
                                  null.model=F)
    
    obj_null<-generalEpis.model(trait.name,
                                scan_pos=m,definScan='parent',
                                cof=sel_cof,definCof=sel_definCof,
                                parGeno,
                                progPheno,
                                megaM,
                                str.residual,
                                null.model=T)
    
    obj_family<-generalEpis.model(trait.name,
                                  scan_pos=m,definScan='family',
                                  cof=sel_cof,definCof=sel_definCof,
                                  parGeno,
                                  progPheno,
                                  megaM,
                                  str.residual,
                                  null.model=F)
    
    
    dev_main_par<--2*(obj_null$asreml.obj$loglik-obj_parent$asreml.obj$loglik)
    (minlog10p_main_par[m]<--log10(0.5*pchisq(dev_main_par,1,lower.tail=FALSE)))
    
    dev_main_fam<--2*(obj_null$asreml.obj$loglik-obj_family$asreml.obj$loglik)
    (minlog10p_main_fam[m]<--log10(0.5*pchisq(dev_main_fam,1,lower.tail=FALSE)))
  
    
    IC.family<- infoCriteria.asreml(obj_family$asreml.obj)
    IC.parent<- infoCriteria.asreml(obj_parent$asreml.obj)
    (AIC_epis[m]<-(IC.parent$AIC-IC.family$AIC))
    
    if (m %% 25 == 0) {
      print(m)
    }
  }
  
  df<- data.frame(ndx=1:nloc,minlog10p_main_fam,minlog10p_main_par,AIC_epis,QTLregion = QTL_region, trait=trait.name)
  df<-cbind(map,df)
  
  minlog10p_main=c()
  
  cF<-which(df$AIC_epis>AIC.threshold) # choose family mod
  cP<-which(!df$AIC_epis>AIC.threshold) # choose parent mod
  
  minlog10p_main[cF]<-df$minlog10p_main_fam[cF]
  minlog10p_main[cP]<-df$minlog10p_main_par[cP]
  
  QTLdefin<-c()
  QTLdefin[cF]<-'family'
  QTLdefin[cP]<-'parent'
  
  df<-df %>% mutate(minlog10p=minlog10p_main,QTLdefin=QTLdefin)
  
  return(df)
}

selEpisQTL.IBD<-function(parGeno,
                         progHap,
                         progPheno,
                         map,
                         trait.name,
                         QTLwindow,
                         threshold=3,
                         str.residual='hete',
                         MQM=TRUE,
                         AIC.threshold=0){
  
  parents = names(parGeno)[-c(1:3)]
  nloc = nrow(map)
  megaM=design.HapM(progHap,parents,1:nloc)
  cofactors=NULL
  definCof=NULL
  coflist=list(cofactors=cofactors,definCof=definCof)
  for (i in 1:100) {
    
    df.scan<-CIMscan.genEpis(trait.name,
                             megaM,
                             parGeno,
                             progPheno,
                             map,
                             cof=coflist[['cofactors']],
                             definCof=coflist[['definCof']],
                             QTLwindow,
                             str.residual,
                             AIC.threshold)
    
    plotQTLscan(df.scan,threshold=threshold,coflist[['cofactors']],trait.name=trait.name)
    
    df.scan.sel = filter(df.scan, QTLregion==FALSE)
    # remove NAs from minlog10p values
    df.scan.sel = filter(df.scan.sel, is.na(minlog10p) == FALSE)
    ord = order(df.scan.sel$minlog10p, decreasing = TRUE)
    max_value = df.scan.sel[ord[1],]$minlog10p
    
    if (max_value < threshold) {
      if (!is.null(cofactors)) cofactors <- sort(cofactors)
    }
    
    if (max_value < threshold) break
    
    if (MQM==F) break
    # if (length(cofactors)>3) break
  
    newcof<-df.scan.sel[ord[1],]$ndx
    cofactors = c(cofactors, newcof)
    definCof=c(definCof,df.scan$QTLdefin[newcof]) 
    coflist<-list(cofactors=cofactors,definCof=definCof)
    
  }
  results <- list(coflist, df.scan)
  return(results)
}

#Calculate goodness-of-fit by random parameters----

## the defination of average semivariance (ASV) see (Piepho 2019, equation 11)
calcu.ASV<-function(V){ 
  n=nrow(V)
  P_lambda=diag(n)-matrix(1/n,n,n)
  ASV_V<-1/(n-1)*psych::tr(V%*%P_lambda)
  return(ASV_V)
}

##the function to calculate goodness-of-fit by random effects (Piepho 2019)
calcR_sq_random<-function(candQTL,trait.name,parGeno,progPheno,progHap,map){
  parents = names(parGeno)[-c(1:3)]
  nP<-length(parents)
  nloc = nrow(map)
  megaM=design.HapM(progHap,parents,1:nloc)
  
  nF<-length(unique(progPheno$pop))
  n.QTL<-length(candQTL$cofactors)
  
#____________________________________________________________________________________  
  ##get the asreml object in the final model fitting QTL candidates and their definations (parent-specific or family-specific)
  if (length(candQTL$cofactors)>1) {
    obj_final<-generalEpis.model(trait.name,
                                 scan_pos=candQTL$cofactors[1],definScan=candQTL$definCof[1],
                                 cof=candQTL$cofactors[-1],definCof=candQTL$definCof[-1],
                                 parGeno,progPheno,megaM,
                                 str.residual='hete',
                                 null.model=F)
  } else {
    obj_final<-generalEpis.model(trait.name[t],
                                 scan_pos=candQTL$cofactors[1],definScan=candQTL$definCof[1],
                                 cof=NULL,definCof=NULL,
                                 parGeno4map_selF,progPheno4map_selF,megaM,
                                 str.residual='hete',
                                 null.model=F)
  }
  
  summary(obj_final$asreml.obj)$varcom # the variance components of QTL candidates
  
  ## construct the G matrix (full model) whose elements are genetic varinces from family- or/and parent- specific effects

  GMat.diag<-c()
  for (q in 1:length(candQTL$cofactors)) {

    if(candQTL$definCof[q]=='parent'){
      varcom.int<-grep(paste0('M',candQTL$cofactors[q],'p'),
                       rownames(summary(obj_final$asreml.obj)$varcom))
      G.component<-summary(obj_final$asreml.obj)$varcom[varcom.int,]$component
      G<-rep(G.component,each= nP)
      GMat.diag<-c(GMat.diag,G)
    } 
    
    if(candQTL$definCof[q]=='family') {
      varcom.int<-grep(paste0('M',candQTL$cofactors[q],'f'),
                       rownames(summary(obj_final$asreml.obj)$varcom))
      G.component<-summary(obj_final$asreml.obj)$varcom[varcom.int,]$component
      G<-rep(G.component,each= 2)
      GMat.diag<-c(GMat.diag,G)
    }
    
  }
  
  (GMat<-as.matrix(diag(GMat.diag)))
  
##construct the Z matrix which is the combination of design matrix from family- or/and parent- specific QTLs
  for (q in 1:length(candQTL$cofactors)) {
    
    if(candQTL$definCof[q]=='parent'){
      z.int<-grep(paste0('M',candQTL$cofactors[q],'p'),names(obj_final$data))
      zmat<-obj_final$data[,z.int]
    }
    
    if(candQTL$definCof[q]=='family'){
      z.int<-grep(paste0('M',candQTL$cofactors[q],'f'),names(obj_final$data))
      zmat<-obj_final$data[,z.int]
    }
    if (q==1){
      allzmat<-zmat
    } else {
      allzmat<-cbind(allzmat,zmat)
    }
  }
  
  allzmat<-as.matrix(allzmat)
  
  ## the ASV of Full model
  V_QTL<-allzmat%*%GMat%*%t(allzmat)
  ASV_QTL<-calcu.ASV(V_QTL)
  
#____________________________________________________________________________________  
  ## get the asreml object from NULL model fitting only population effects when set null.model=T
  
  obj_null<-generalEpis.model(trait.name,
                              scan_pos=candQTL$cofactors[1],definScan=candQTL$definCof[1],
                              cof=NULL,definCof=NULL,
                              parGeno,progPheno,megaM,
                              str.residual='hete',
                              null.model=T)
  
  
  summary(obj_null$asreml.obj)$varcom
  R.component<-summary(obj_null$asreml.obj)$varcom$component
  
  #family-specific residual matrix
  V_0<-as.matrix(diag(rep(R.component,as.numeric(table(progPheno$pop)))))
  ASV_0<-calcu.ASV(V_0)
  
  ## R^2 of Random effect
  R_sq_random=ASV_QTL/ASV_0
  
  return(R_sq_random)
}



