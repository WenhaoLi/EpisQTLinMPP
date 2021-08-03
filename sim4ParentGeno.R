nmar.chr<- c(100,100,100,100,100)
n.chr<-length(nmar.chr)
len.chr<- c(100,100,100,100,100)
npar<- 4

name.mar<-paste0('M',1:sum(nmar.chr))

allname.chr=c()
for (c in 1:n.chr) {
  name.chr<-rep(c(1:n.chr)[c],each=nmar.chr[c])
  allname.chr<-c(allname.chr,name.chr)
}

allpos.chr=c()
for (c in 1:n.chr) {
  pos.chr<-seq(0.1,len.chr[c],length.out = nmar.chr[c])
  allpos.chr<-c(allpos.chr,pos.chr)
}

allpos.chr<-round(allpos.chr,1)
parmarker<-data.frame(marker=name.mar,chromosome=allname.chr,pos=allpos.chr)

pargeno<-as.data.frame(matrix(sample(c(1,2),size=sum(nmar.chr)*npar,replace = T),ncol = npar))
names(pargeno)<-paste0('parent_',LETTERS[1:npar])

pargenotype<-cbind(parmarker,pargeno)
# write.csv(pargenotype,'pargenotype5chr.csv',row.names = F)
