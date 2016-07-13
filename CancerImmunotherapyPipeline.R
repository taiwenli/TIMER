## License GPL2.0
## Developed by Bo Li, bli@jimmy.harvard.edu, 2016
## TIMER Pipeline for analyzing immune cell components in the tumor microenvironment

## Before running this code, make sure you have all the datasets available in the right format:
##  1. Gene expression data for all 23 tumors, could be directly downloaded from GDAC firehose, but needs to be in .Rdata format.
##  2. Clinical data for all samples, directly downloaded from TCGA public ftp
##  3. Auxilliary datasets available from TIMER website


##----- setup parameters and establish the output file -----##
if(cc=='skcm')cc.type='06A' else cc.type='01A'
cur.dir='/Users/bo/Desktop/DFCI/microEnv/'  ## Change this to your own analysis directory
setwd(paste(cur.dir,'results',cc,sep='/'))
options(scipen=6)
write(paste(cc,' output\n'),file=paste(cur.dir,'/results/',cc,'/output-statistics.txt',sep=''))

signature.genes=c('CD19','TRAT1','CD8B','CCR3','CD163','CCL17')
names(signature.genes)=c('B_cell','T_cell.CD4','T_cell.CD8','Neutrophil','Macrophage','DC')

##----- load and process gene expression data -----##
if(cc %in% c('gbm','ov')){
  dd=get(load(paste(cur.dir,'data/RNAseq/',cc,'Affy.Rdata',sep='')))
}else dd=get(load(paste(cur.dir,'data/RNAseq/',cc,'RNAseq.Rdata',sep='')))
dd=as.matrix(dd)
mode(dd)='numeric'
if(!cc %in% c('gbm','ov','esca','stad'))dd=dd*1e6   ## rsem scaled estimates needs multiply 1e6, Array or RPKM does not need.
tmp=strsplit(rownames(dd),'\\|')
tmp=sapply(tmp,function(x)x[[1]])
tmp.vv=which(nchar(tmp)>1)
rownames(dd)=tmp
dd=dd[tmp.vv,]

##----- load immune marker genes from Abbas et al., 2005 -----##
tmp=read.csv(paste(cur.dir,'data/immune_datasets/IRIS-marker-gene.txt',sep=''),header=T,sep='\t',stringsAsFactors=F)
marker.list=tmp[,1]
names(marker.list)=tmp[,7]
names(marker.list)=gsub(' ','_',tmp[,7])
names(marker.list)=gsub('Dendritic_Cell','DC',names(marker.list))
names(marker.list)=gsub('Neutrophil','Neutrophils',names(marker.list))
gsub('Cell','cell',names(marker.list))->names(marker.list)

##----- load reference data of sorted immune cells -----##
load(paste(cur.dir,'data/immune_datasets/HPCTimmune.Rdata',sep=''))

##----- load and process tumor purity data -----##
AGP=read.table(paste(cur.dir,'data/AGP/AGP-',cc,'.txt',sep=''),sep='\t',header=T)
AGP=AGP[which(AGP[,'PoP']>0.01),]
tmp=strsplit(rownames(HPCT.immune),';')
AffyIDtoGenes=sapply(tmp,function(x)x[[1]])
names(AffyIDtoGenes)=sapply(tmp,function(x)x[[2]])
marker.list.genes=AffyIDtoGenes[marker.list]

##----- function to edit TCGA ID, with the option of keeping the first num.res fields -----##
getID <- function(sID,num.res=3){
  mm=c()
  for(id in sID){
    tmp=unlist(strsplit(id,'-'))
    if(length(tmp)==1){
      tmp=unlist(strsplit(id,'\\.'))
    }
    ll='TCGA'
    for(j in 2:num.res){
      ll=paste(ll,tmp[j],sep='-')
    }
    mm=c(mm,ll)
  }
  return(mm)
}
rownames(AGP)=getID(AGP[,1],4)
colnames(dd)=getID(colnames(dd),4)

##----- Select single reference samples of pre-selected immune cell types -----##
B_cell=362:385
T_cell.CD4=grep('T_cell.CD4',colnames(HPCT.immune))
T_cell.CD8=grep('T_cell.CD8',colnames(HPCT.immune))
NK=328:331
Neutrophil=344:361
Macrophage=66:80
DC=151:238
curated.ref=HPCT.immune[,c(B_cell,T_cell.CD4,T_cell.CD8,NK,Neutrophil,Macrophage,DC)]

curated.cell.types=colnames(curated.ref)
names(curated.cell.types)=c(rep('B_cell',length(B_cell)),rep('T_cell.CD4',length(T_cell.CD4)),rep('T_cell.CD8',length(T_cell.CD8)),rep('NK',length(NK)),rep('Neutrophil',length(Neutrophil)),rep('Macrophage',length(Macrophage)),rep('DC',length(DC)))

##----- function to preprocess the reference dataset, not necessary if the processed data "curated.ref.genes.Rdata" is available -----##
#GetRefGenes <- function(curated.ref){
#  AffyIDwGenes= rownames(HPCT.immune)
#  library(sqldf)
#  tmpDD=data.frame(curated.ref)
#  tmpDD=tmpDD[order(rownames(tmpDD)),]
#  colnames(tmpDD)=gsub('\\.','_',colnames(tmpDD))
#  genes=sapply(strsplit(rownames(tmpDD),';'),function(x)x[[1]])
#  tmpDD=cbind(genes,tmpDD)
#  tmpDD=tmpDD[order(genes),]
#  tmp0=c()
#  for(i in colnames(tmpDD)[2:ncol(tmpDD)]){
#    print(i)
#    tmp=sqldf(paste('select max(',i,') from tmpDD group by genes',sep=''))
#    if(length(tmp0)==0)tmp0=tmp else tmp0=cbind(tmp0,tmp)
#  }
#  colnames(tmp0)=colnames(tmpDD)[2:ncol(tmpDD)]
#  rownames(tmp0)=unique(tmpDD[,1])
#  curated.ref.genes=tmp0
#  return(curated.ref.genes)
#}

load(paste(cur.dir,'data/immune_datasets/curated.ref.genes.Rdata',sep=''))

##----- Combine TCGA gene expression profiles with the selected reference data, remove batch effect and aggregate samples of each immune category by taking the median -----##
RemoveBatchEffect <- function(){
  library(sva)
  tmp.dd=as.matrix(dd)
  tmp=sapply(strsplit(rownames(dd),'\\|'),function(x)x[[1]])
  rownames(tmp.dd)=tmp
  tmp.dd=tmp.dd[which(nchar(tmp)>1),]
  tmp.ss=intersect(rownames(tmp.dd),rownames(curated.ref.genes))
  N1=ncol(tmp.dd)
  tmp.dd=cbind(tmp.dd[tmp.ss,],curated.ref.genes[tmp.ss,])
  tmp.dd=as.matrix(tmp.dd)
  mode(tmp.dd)='numeric'
  N2=ncol(curated.ref.genes)
  tmp.batch=c(rep(1,N1),rep(2,N2))
  tmp.dd0=ComBat(tmp.dd,tmp.batch,c())
  dd.br=tmp.dd0[,1:N1]
  curated.ref.genes.br=tmp.dd0[,(N1+1):(N1+N2)]
  tmp0=c()
  for(kk in unique(names(curated.cell.types))){
    tmp.vv=which(names(curated.cell.types)==kk)
    tmp0=cbind(tmp0,apply(curated.ref.genes.br[,tmp.vv],1,median,na.rm=T))
  }
  curated.ref.genes.agg.br=tmp0
  colnames(curated.ref.genes.agg.br)=unique(names(curated.cell.types))
  #rownames(curated.ref.genes.agg.br)=rownames(curated.ref.genes.br)
  return(list(dd=dd.br,rr=curated.ref.genes.br,rrg=curated.ref.genes.agg.br))
}

tmp=RemoveBatchEffect()
dd.br=tmp$dd
curated.ref.genes.br=tmp$rr
curated.ref.genes.agg.br=tmp$rrg

##----- Constrained regression method implemented in Abbas et al., 2009 -----##
getFractions.Abbas <- function(XX,YY,w=NA){
  ss.remove=c()
  ss.names=colnames(XX)
  while(T){
    if(length(ss.remove)==0)tmp.XX=XX else{
      if(is.null(ncol(tmp.XX)))return(rep(0,ncol(XX)))
      tmp.XX=tmp.XX[,-ss.remove]
    }
    if(length(ss.remove)>0){
      ss.names=ss.names[-ss.remove]
      if(length(ss.names)==0)return(rep(0,ncol(XX)))
    }
    if(is.na(w[1]))tmp=lsfit(tmp.XX,YY,intercept=F) else tmp=lsfit(tmp.XX,YY,w,intercept=F)
    if(is.null(ncol(tmp.XX)))tmp.beta=tmp$coefficients[1] else tmp.beta=tmp$coefficients[1:(ncol(tmp.XX)+0)]
    if(min(tmp.beta>0))break
    ss.remove=which.min(tmp.beta)
  }
  tmp.F=rep(0,ncol(XX))
  names(tmp.F)=colnames(XX)
  tmp.F[ss.names]=tmp.beta
  return(tmp.F)
}

##----- function to calculate the residuals from regression -----##
fn <- function(beta0,XX,Y)return(log(sum(abs(Y-XX%*%beta0))))

##----- function to select genes with expression values negatively correlated with tumor purity -----##
getPurityGenes <- function(dd,AGP,thr.p=0.05,thr.c=0,mode='env'){
  tmp.ss=intersect(colnames(dd),rownames(AGP))
  if(length(tmp.ss)==0){
    colnames(dd)=getID(colnames(dd))
    tmp.ss=intersect(colnames(dd),rownames(AGP))
  }
  tmp.dd=dd[,tmp.ss]
  tmp=lapply(rownames(tmp.dd),function(x)cor.test(tmp.dd[x,],as.numeric(AGP[colnames(tmp.dd),2]),method='s'))
  tmp.pp=sapply(tmp,function(x)x$p.value)
  tmp.cor=sapply(tmp,function(x)x$estimate)
  names(tmp.pp)=names(tmp.cor)=rownames(dd)
  if(mode=='env')vv=names(which(tmp.pp <=thr.p&tmp.cor < thr.c))
  if(mode=='tumor')vv=names(which(tmp.pp <=thr.p&tmp.cor > thr.c))
  return(vv)
}

##----- selection genes negatively correlated with purity and overlap with immune marker genes -----##
vv.t=getPurityGenes(dd,AGP,thr.p=0.05,thr.c= -0.2)
vv.t=intersect(vv.t,rownames(curated.ref.genes.agg.br))
vv=intersect(vv.t,marker.list.genes)

##----- remove outlier genes whose expression may drive the colinearity of similar covariates in the regression -----##
RemoveOutliers <- function(vv, ref.dd, thr.q=0.99){
  ## removes upper thr.q quantile for every reference feature
  remove.vv=c()
  for(i in 1:ncol(ref.dd)){
    tmp=quantile(ref.dd[vv,i],thr.q)[1]
    tmp.vv=which(ref.dd[vv,i]>tmp)
    remove.vv=c(remove.vv,tmp.vv)
  }
  remove.vv=unique(remove.vv)
  return(vv[-remove.vv])
}

##---- calculate differences between the correlations of reference immune cells using Pearson's or Spearman's correlations -----##
tmp.diff=sum(sum(abs(cor(curated.ref.genes.agg.br[vv,],method='p')-cor(curated.ref.genes.agg.br[vv,],method='s'))))

if(tmp.diff>= -10000){
  vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,-4])
  vv=vv0
}

cat("Number of genes inversely correlated with purity is ",length(vv.t),'\n\n',sep='',file='output-statistics.txt',append=T)
cat("Number of immune genes inversely correlated with purity is ",length(vv),'\n\n',sep='',file='output-statistics.txt',append=T)

##----- calculate the significance of enrichment for purity selected genes to immune marker genes -----##
tmp.ss0=intersect(rownames(curated.ref.genes.agg.br),rownames(dd.br))
n.immune=length(intersect(marker.list.genes,tmp.ss0))
cat("Test if immune genes are enriched for inverse correlation with purity: \n\n",file='output-statistics.txt',append=T)
sink(file='output-statistics.txt',append=T);print(fisher.test(matrix(c(length(vv),length(vv.t)-length(vv),n.immune,length(tmp.ss0)-n.immune),2,2)));sink()

##----- function to process deconvolution method in batch -----##
BatchFractions <- function(XX,YYd){
  Fmat=c()
  for(i in 1:ncol(YYd)){
    YY=YYd[,i]
    tmp.F=getFractions.Abbas(XX,YY)
    #tmp.F=getFractions.Optim(XX,YY)
    Fmat=rbind(Fmat,tmp.F)
  }
  rownames(Fmat)=colnames(YYd)
  colnames(Fmat)=colnames(XX)
  return(Fmat)
}

##----- perform batch deconvolution -----##
XX=curated.ref.genes.agg.br[vv,c(-4)]
YYd=dd.br[vv,]
Fmat=BatchFractions(XX,YYd)

##----- CD4 and CD8 T cells are likely to be similar, resulting in colinearity. Codes below are procedures to remove outlier genes that may result in colinearity until the two covariates are linearly separable. -----##
if(cor(Fmat[,2],Fmat[,3])<= -0.2){
  if(tmp.diff>=1){
    tmp.cor=c()
    thr.qlist=c(0.99)
    for(tq in thr.qlist){
      vv=intersect(vv.t,marker.list.genes)
      vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,-4],tq)
      vv=vv0
      XX=curated.ref.genes.agg.br[vv,]
      YYd=dd.br[vv,]
      tmp.Fmat=BatchFractions(XX,YYd)
      tmp.cor=c(tmp.cor,cor(tmp.Fmat[,2],tmp.Fmat[,3],method='s'))
    }
    tmp.vv=which.max(tmp.cor)
    tq=thr.qlist[tmp.vv]
    vv=intersect(vv.t,marker.list.genes)
    vv0=RemoveOutliers(vv,curated.ref.genes.agg.br[,c(-4)],tq)
    vv=vv0
    XX=curated.ref.genes.agg.br[vv,c(-4)]
    YYd=dd.br[vv,]
    Fmat=BatchFractions(XX,YYd)
    Fmat0.p=Fmat[grep(cc.type,rownames(Fmat)),]
    rownames(Fmat0.p)=getID(rownames(Fmat0.p))
  }
}

while(cor(Fmat[,2],Fmat[,3])<=-0.3){
  if(length(vv)<=50)break
  vv=vv[-as.numeric(names(table(apply(dd[vv,],2,which.max))))]
  XX=curated.ref.genes.agg.br[vv,c(-4)]
  #XX=XX[,!colnames(XX) %in% tmp.remove]
  YYd=dd.br[vv,]
  Fmat=BatchFractions(XX,YYd)
}

Fmat0.p=Fmat[grep(cc.type,rownames(Fmat)),]
rownames(Fmat0.p)=getID(rownames(Fmat0.p))


##----- Simulation based evaluation of inference based on 1) gene set vv and 2) gene-gene correlation structure of the cancer data -----##
library(MASS)
thr.sd=1
N=nrow(Fmat)
N=1
XX=curated.ref.genes.agg.br[vv,-4]
sim.Fmat=c()
est.Fmat=c()
sim.YY=c()
SIGMA=cov(t(dd.br[vv,]))
for(i in 1:N){
    tmp=runif(6)
    tmp=tmp/sum(tmp)
    mm=XX %*% tmp
    mm=mm+rnorm(length(mm),0,thr.sd)
    yy=mvrnorm(1,mm,SIGMA)
    sim.Fmat=rbind(sim.Fmat,tmp)
    sim.YY=cbind(sim.YY,yy)
}
rownames(sim.YY)=vv
est.Fmat=BatchFractions(XX,sim.YY)

sim.COR=diag(cor(est.Fmat,sim.Fmat,method='s'))
names(sim.COR)=colnames(XX)
sim.COR=list(sim.COR)
names(sim.COR)=cc

##----- Cytolityc score -----##
CYT=list(log(apply(dd[c('GZMA','PRF1'),],2,mean)))
names(CYT)=cc

####----------------------------------------#####
## Below are analysis part: run in tandem with above codes.

## Quality Control: All should satisfy. If any is violated, the results cannot be trusted.
#1. Check if distribution of reference dataset can be linearly unified with tumor expression
pdf(paste(cur.dir,'/results/',cc,'/output.pdf',sep=''))
par(cex.lab=1,cex.axis=1,las=2,pch=19,mar=c(5,5,2,0))
qqplot(dd.br[,1],XX[,1],xlab='Tumor Expression',ylab='Ref Expression',main='Sample-Sample Q-Q plot') ## q-q plot by sample should look like a straight line. Extreme values may saturate for Affy array data, but most of the data should align well.
qqplot(dd.br[1,],curated.ref.genes[1,],xlab='Tumor Expression',ylab='Ref Expression',main='Gene-Gene Q-Q plot') ## q-q plot by gene

#2. Check co-linearity of immune cell expression
plot(data.frame(XX),pch=19,cex=0.2) ## reference colinearity
plot(data.frame(Fmat),pch=19,cex=0.2,col=2) ## estimation correlations

#3. Check correlations for predicted fractions
load('/Users/bo/Desktop/bo/redblue.Rdata')
heatmap(cor(Fmat),col=redblue)

hist(cor(Fmat)[upper.tri(cor(Fmat))]) ## negative correlations are driven by co-linearity. Positive correlations are likely to be true.

#4. Check if fractions of individual immune cell types are correlated with the marker gene expression for that cell type
par(mfrow=c(2,2))
cat('QC4: Correlation of Estimated Abundance with Marker Genes\n\n',file='output-statistics.txt',append=T)

#B cell
if('B_cell' %in% colnames(Fmat)){
  plot(log(dd.br['CD19',]),Fmat[,'B_cell'],xlab='CD19 expression',ylab='B-cell abundance',pch=19)
  cat('CD19\n\n',file='output-statistics.txt',append=T)
  sink(file='output-statistics.txt',append=T);print(cor(log(dd.br['CD19',]),Fmat,use='complete.obs'));sink()
}
#T CD8B
if('T_cell.CD8' %in% colnames(Fmat)){
  plot(log(dd.br['CD8B',]),Fmat[,'T_cell.CD8'],xlab='CD8B expression',ylab='T-cell.CD8 abundance',pch=19)
  cat('CD8B\n\n',file='output-statistics.txt',append=T)
  sink(file='output-statistics.txt',append=T);print(cor(log(dd.br['CD8B',]),Fmat,use='complete.obs'));sink()
}
if('Macrophage' %in% colnames(Fmat)){
  #Macrophage
  plot(log(dd.br['CD163',]),Fmat[,'Macrophage'],xlab='CD163 expression',ylab='Macrophage abundance',pch=19)
  cat('CD163\n\n',file='output-statistics.txt',append=T)
  sink(file='output-statistics.txt',append=T);print(cor(log(dd.br['CD163',]),Fmat,use='complete.obs'));sink()
}
#DC
if('DC' %in% colnames(Fmat)){
  plot(log(dd.br['ITGAM',]),Fmat[,'DC'],xlab='CD11b expression',ylab='Dendritic Cell abundance',pch=19)
  cat('ITGAM\n\n',file='output-statistics.txt',append=T)
  sink(file='output-statistics.txt',append=T);print(cor(log(dd.br['ITGAM',]),Fmat,use='complete.obs'));sink()
}
if('Neutrophil' %in% colnames(Fmat)){
  #Neutrophil
  plot(log(dd.br['AMPD2',]),Fmat[,'Neutrophil'],xlab='AMPD2 expression',ylab='Neutrophil abundance',pch=19)
  cat('AMPD2\n\n',file='output-statistics.txt',append=T)
  sink(file='output-statistics.txt',append=T);print(cor(log(dd.br['AMPD2',]),Fmat,use='complete.obs'));sink()
}
#5. Check if predicted neuthrophil abundance is positively correlated with necrosis
necrosisData=read.csv(paste(cur.dir,'/data/Clinical/nationwidechildrens.org_biospecimen_slide_',cc,'.txt',sep=''),sep='\t',header=T)
necrosisData=as.matrix(necrosisData)
necrosisData=necrosisData[grep(cc.type,necrosisData[,1]),]
rownames(necrosisData)=getID(necrosisData[,1])
tmp.ss=intersect(rownames(necrosisData),rownames(Fmat0.p))
if('Neutrophil' %in% colnames(Fmat)){
  plot(necrosisData[tmp.ss,8],Fmat0.p[tmp.ss,'Neutrophil'],xlab='%Necrosis',ylab='Neutrophil Abundance')
  cat('\n\n-------necrosis vs neutrophil-------\n\n',file='output-statistics.txt',append=T)
  sink(file='output-statistics.txt',append=T);print(cor.test(as.numeric(necrosisData[tmp.ss,8]),Fmat0.p[tmp.ss,'Neutrophil'],use='complete.obs'));sink()
}

## Survival Analysis
library(survival)
dd.surv=read.csv(paste(cur.dir,'data/Clinical/nationwidechildrens.org_clinical_patient_',cc,'.txt',sep=''),header=T,sep='\t',stringsAsFactors=F)
dd.surv=dd.surv[3:nrow(dd.surv),]
rownames(dd.surv)=dd.surv[,1]
dat.surv=matrix(NA,ncol=2,nrow=nrow(dd.surv))
tmp=dd.surv[,'vital_status']
tmp=as.numeric(as.factor(tmp))-1
dat.surv[,1]=tmp
mode(dat.surv)='numeric'
rownames(dat.surv)=rownames(dd.surv)
tmp.vv0=grep('last_followup',colnames(dd.surv))
if(length(tmp.vv0)==0)tmp.vv0=grep('last_contact',colnames(dd.surv))
dat.surv[which(tmp==0),2]=as.numeric(dd.surv[which(tmp==0),tmp.vv0])
tmp.vv1=grep('death',colnames(dd.surv))
if(length(tmp.vv1)>1){
  tmp.vv=grep('days',colnames(dd.surv))
  tmp.vv1=intersect(tmp.vv1,tmp.vv)
}
dat.surv[which(tmp==1),2]=as.numeric(dd.surv[which(tmp==1),tmp.vv1])
dd.surv=as.matrix(dd.surv)
tmp.vv=intersect(grep('age_',colnames(dd.surv)),grep('diagnosis',colnames(dd.surv)))
cc.ages=as.numeric(dd.surv[,tmp.vv])
names(cc.ages)=rownames(dd.surv)
tmp.vv=c(grep('tumor_stage',colnames(dd.surv)),grep('clinical_stage',colnames(dd.surv)))
if(length(tmp.vv)>0){
  if(length(tmp.vv)>1)tmp.vv=tmp.vv[1]
  cc.stages=dd.surv[,tmp.vv]
  names(cc.stages)=rownames(dd.surv)
  cc.stages[which(cc.stages%in%c('Stage I','Stage IA','Stage IB'))]='Stage-I'
  cc.stages[which(cc.stages%in%c('Stage II','Stage IIA','Stage IIB'))]='Stage-II'
  cc.stages[-grep('-',cc.stages)]='Stage-II+'
  cc.stages[which(cc.stages %in% names(which(table(cc.stages)<=15)))]=NA
}else cc.stages=c()
if(length(unique(na.omit(cc.stages)))==1)cc.stages=c()
tmp.vv=grep('grade',colnames(dd.surv))
if(length(tmp.vv)>0){
  cc.grades=dd.surv[,tmp.vv[1]]
  names(cc.grades)=rownames(dd.surv)
  cc.grades[which(cc.grades=='GX'|cc.grades=='[Not Available]')]=NA
  cc.grades[which(cc.grades %in% names(which(table(cc.grades)<=15)))]=NA
  if(length(unique(na.omit(cc.grades)))==1)cc.grades=c()
}else cc.grades=c()

tmp.ss=intersect(rownames(dat.surv),rownames(Fmat0.p))
dat.surv[which(dat.surv[,1]==2),1]=NA
tmp.Surv=Surv(dat.surv[tmp.ss,2],dat.surv[tmp.ss,1])

cat('\n\n------Correlation with patient survival--------\n\n',file='output-statistics.txt',append=T)
if(length(cc.stages)>0){
  for(i in 1:ncol(Fmat0.p)){
    cat('\n',colnames(Fmat0.p)[i],'\n\n',file='output-statistics.txt',append=T)
    sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat0.p[tmp.ss,i]+cc.ages[tmp.ss]+cc.stages[tmp.ss]));sink()
  }
  if('DC' %in% colnames(Fmat)){
    sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat0.p[tmp.ss,'Macrophage']+Fmat0.p[tmp.ss,'DC']+cc.ages[tmp.ss]+cc.stages[tmp.ss]));sink()
  }
  sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat0.p[tmp.ss,]+cc.ages[tmp.ss]+cc.stages[tmp.ss]));sink()
}else{
  for(i in 1:ncol(Fmat0.p)){
    cat('\n',colnames(Fmat0.p)[i],'\n\n',file='output-statistics.txt',append=T)
    sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat0.p[tmp.ss,i]+cc.ages[tmp.ss]));sink()
  }
  if('DC' %in% colnames(Fmat)){
    sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat0.p[tmp.ss,'Macrophage']+Fmat0.p[tmp.ss,'DC']+cc.ages[tmp.ss]));sink()
  }
  sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat0.p[tmp.ss,]+cc.ages[tmp.ss]));sink()
}

## smoking
tmp.vv=grep('tobacco_smoking_history_indicator',colnames(dd.surv))
if(length(tmp.vv)>0){
  smoking=dd.surv[,tmp.vv]
  smoking[grep('\\[',smoking)]=NA
  smoking[grep('Not Specified',smoking)]=NA
  names(smoking)=rownames(dd.surv)
  tmp=factor(smoking,levels=c('Current smoker','Current reformed smoker for < or = 15 years','Current reformed smoker for > 15 years','Lifelong Non-smoker'))
  smoking=tmp
  names(smoking)=rownames(dd.surv)
  cat("\n\n-------Correlation with smoking--------\n\n",file='output-statistics.txt',append=T)
  if(!is.null(cc.stages)){sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~smoking[tmp.ss]+cc.ages[tmp.ss]+cc.stages[tmp.ss]));sink()}else{sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~smoking[tmp.ss]+cc.ages[tmp.ss]));sink()}
}

cat("\n\n-------Correlation with patient survival: residual after purity regression--------\n\n",file='output-statistics.txt',append=T)

tmp.ss01A=grep(cc.type,rownames(Fmat))
rownames(Fmat)=getID(rownames(Fmat),4)
Fmat.res=c()
for(i in 1:ncol(Fmat)){
  tmp=intersect(rownames(Fmat[tmp.ss01A,]),rownames(AGP))
  tmp.res=lm(Fmat[tmp,i]~AGP[tmp,2])$residuals
  Fmat.res=cbind(Fmat.res,tmp.res)
}

cat("\n\n............Find Genes highly expressed in tumor correlated with immune cell abundance............\n\n",file='output-statistics.txt',append=T)

colnames(dd.br)=getID(colnames(dd.br),4)
tumor.vv=getPurityGenes(dd,AGP,thr.p=0.05,thr.c= 0.2,mode='tumor')
tumor.vv=intersect(tumor.vv,rownames(dd.br))
dd.res=matrix(NA,nrow=length(tumor.vv),ncol=nrow(Fmat.res))
rownames(dd.res)=tumor.vv
colnames(dd.res)=rownames(Fmat.res)
colnames(Fmat.res)=colnames(Fmat)
for(k in tumor.vv){
  tmp=intersect(rownames(Fmat[tmp.ss01A,]),rownames(AGP))
  tmp.res=lm(dd.br[k,tmp]~AGP[tmp,2])$residuals
  dd.res[k,names(tmp.res)]=tmp.res
}

#tumor.vv=getPurityGenes(dd,AGP,thr.p=0.05,thr.c= 0.3,mode='tumor')
#tumor.vv=intersect(tumor.vv,rownames(dd.br))
top.gg=c()
tmp.names=c()
for(i in 1:ncol(Fmat)){
  tmp.sp=intersect(rownames(Fmat),rownames(AGP))
  if(cor(Fmat[tmp.sp,i],AGP[tmp.sp,2],use='complete.obs')> -0.1)next
  tmp.names=c(tmp.names,colnames(Fmat)[i])
  cat('Top correlated tumor genes for',colnames(Fmat.res)[i],'\n\n',file='output-statistics.txt',append=T)
  temp=lapply(1:nrow(dd.res),function(x)cor.test(dd.res[x,],Fmat.res[,i],use='complete.obs'))
  if(length(tumor.vv)>1)tmp.cor0=cor(Fmat[,i],t(dd.br[tumor.vv,])) else tmp.cor0=cor(Fmat[,i],dd.br[tumor.vv,])
  tmp.cor=sapply(temp,function(x)x$estimate)
  tmp.p=sapply(temp,function(x)x$p.value)
  tmp.p=p.adjust(tmp.p)
  names(temp)=names(tmp.cor)=tumor.vv
  top.gg=c(top.gg,list(tmp.cor[which(tmp.p<=0.05&tmp.cor>0 & tmp.cor0<=0.1)]))
}
names(top.gg)=tmp.names

rownames(Fmat.res)=getID(names(tmp.res))
colnames(Fmat.res)=colnames(Fmat)
tmp.ss=intersect(rownames(Fmat.res),rownames(dd.surv))
tmp.Surv=Surv(dat.surv[tmp.ss,2],dat.surv[tmp.ss,1])

## Stage

if(length(cc.stages)>0){
  for(i in 1:ncol(Fmat.res)){
    cat('\n',colnames(Fmat.res)[i],'\n\n',file='output-statistics.txt',append=T)
    sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat.res[tmp.ss,i]+cc.ages[tmp.ss]+cc.stages[tmp.ss]));sink()
  }
  if('DC' %in% colnames(Fmat)){
    sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat.res[tmp.ss,'Macrophage']+Fmat.res[tmp.ss,'DC']+cc.ages[tmp.ss]+cc.stages[tmp.ss]));sink()
  }
  sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat.res[tmp.ss,]+cc.ages[tmp.ss]+cc.stages[tmp.ss]));sink()
}else{
  for(i in 1:ncol(Fmat.res)){
    cat('\n',colnames(Fmat.res)[i],'\n\n',file='output-statistics.txt',append=T)
    sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat.res[tmp.ss,i]+cc.ages[tmp.ss]));sink()
  }
  if('DC' %in% colnames(Fmat)){
    sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat.res[tmp.ss,'Macrophage']+Fmat.res[tmp.ss,'DC']+cc.ages[tmp.ss]));sink()
  }
  sink(file='output-statistics.txt',append=T);print(coxph(tmp.Surv~Fmat.res[tmp.ss,]+cc.ages[tmp.ss]));sink()
}

cat('\n\n--------Correlation with patient age--------\n\n',file='output-statistics.txt',append=T)
sink(file='output-statistics.txt',append=T);print(cor(Fmat0.p[tmp.ss,],cc.ages[tmp.ss],use='complete.obs'));sink()

cat('\n\n--------Correlation with clinical stage--------\n\n',file='output-statistics.txt',append=T)
if(length(cc.stages)>0){
  for(i in 1:ncol(Fmat0.p)){
    cat('\n',colnames(Fmat0.p)[i],'\n\n',file='output-statistics.txt',append=T)
    sink(file='output-statistics.txt',append=T);print(summary(aov(lm(Fmat0.p[tmp.ss,i]~as.factor(cc.stages[tmp.ss])))));sink()
  }
}
## optional analysis on metastasis
# recurrent or metastasis
tmp.vv=grep('tumor_status',colnames(dd.surv))
if(length(tmp.vv)>0){
  par(mfrow=c(2,3))
  cat('\n\n-------tumor_status------\n\n',file='output-statistics.txt',append=T)
  tumor_status=dd.surv[,tmp.vv]
  names(tumor_status)=rownames(dd.surv)
  tumor_status[grep('\\[',tumor_status)]=NA
  for(i in 1:ncol(Fmat0.p)){
    boxplot(Fmat0.p[tmp.ss,i]~tumor_status[tmp.ss])
    t.test(Fmat0.p[tmp.ss,i]~tumor_status[tmp.ss])
    tumor_status.c=as.numeric(as.factor(tumor_status))
    names(tumor_status.c)=names(tumor_status)
    cat('\n',colnames(Fmat0.p)[i],'\n\n',file='output-statistics.txt',append=T)
    sink(file='output-statistics.txt',append=T);print(summary(glm(tumor_status.c[tmp.ss]~Fmat0.p[tmp.ss,i]+cc.ages[tmp.ss])));sink()
    sink(file='output-statistics.txt',append=T);print(t.test(Fmat0.p[tmp.ss,i]~tumor_status.c[tmp.ss]));sink()
  }
}
# race
race=c()
tmp.vv=grep('race',colnames(dd.surv))
if(length(tmp.vv)>0 &length(unique(dd.surv[,tmp.vv]))>1){
  cat('\n\n------race-------\n\n',file='output-statistics.txt',append=T)
  par(mfrow=c(2,3))
  if(length(tmp.vv)>1)tmp.vv=tmp.vv[2]
  race=dd.surv[,tmp.vv]
  race[grep('\\[',race)]=NA
  if(length(unique(na.omit(race)))>1){
    names(race)=rownames(dd.surv)
    race[which(race %in% names(which(table(race)<=10)))]=NA
    if(length(unique(na.omit(race)))>1){
      for(i in 1:ncol(Fmat0.p)){
        boxplot(Fmat0.p[tmp.ss,i]~race[tmp.ss])
        cat('\n',colnames(Fmat0.p)[i],'\n\n',file='output-statistics.txt',append=T)
        sink(file='output-statistics.txt',append=T);print(summary(aov(lm(Fmat0.p[tmp.ss,i]~as.factor(race[tmp.ss])))));sink()
      }
    }
  }else race=c()
}

# gender
gender=c()
tmp.vv=grep('gender',colnames(dd.surv))
if(length(tmp.vv)>1)tmp.vv=tmp.vv[1]
gender=dd.surv[,tmp.vv]
gender[grep('\\[',gender)]=NA
names(gender)=rownames(dd.surv)


## Compare immune cell abundance between primary tumor and adjacent tissue
tmp.11A=grep('11A',rownames(Fmat))
tmp.01A=grep(cc.type,rownames(Fmat))
rownames(Fmat)=getID(rownames(Fmat),4)
tmp.types=sapply(strsplit(rownames(Fmat),'-'),function(x)x[[4]])
tmp.tissue=rep(0,length(c(tmp.01A,tmp.11A)))
names(tmp.tissue)=colnames(dd.br)[c(tmp.01A,tmp.11A)]
tmp.tissue[colnames(dd.br)[tmp.11A]]=1
if(length(tmp.11A)>10){
  par(mfrow=c(2,3))
  for(i in 1:ncol(Fmat))boxplot(Fmat[,i]~tmp.types)
  
  cat("\n\n----------primary tumor compared with adjacent tissue-------\n\n",file='output-statistics.txt',append=T)
  for(i in 1:ncol(Fmat)){
    cat('\n',colnames(Fmat)[i],'\n\n',file='output-statistics.txt',append=T)
    sink(file='output-statistics.txt',append=T);print(t.test(Fmat[tmp.01A,i],Fmat[tmp.11A,i]));sink()
  }
}
dev.off()


setwd(paste(cur.dir,'scripts',sep='/'))





























