
BinMeth<-function(data_in){

Input<-as.numeric(data_in)

ProcClust<-kmeans(Input[which(is.na(Input)==F)],2)
Proc<-ProcClust$cluster
low<-which(ProcClust$centers==min(ProcClust$centers))
high<-which(ProcClust$centers==max(ProcClust$centers))
Proc[which(Proc==low)]<-0
Proc[which(Proc==high)]<-1

Output<-Input
Output[which(is.na(Output)==F)]<-Proc
return(Output)
}


CpGall<-read.table("CpG_genes_percentMeth.bed", sep="\t", stringsAsFactors=F,header=T)

cols_to_convert <- names(CpGall)[5:ncol(CpGall)]

CpGall <- CpGall %>%
  mutate_at(vars(cols_to_convert), as.numeric)

BinTab<-CpGall

for(i in 5:ncol(CpGall)){
BinTab[,i]<-BinMeth(CpGall[,i])
}

FCFCmean<-apply(CpGall[,5:7],1,mean,na.rm=T)
NCNCmean<-apply(CpGall[,15:18],1,mean,na.rm=T)

MeanTab<-cbind(CpGall[,1:4],FCFCmean,NCNCmean)






FCMeth<-data.frame()
for(i in 5:7){

a<-which(is.na(BinTab[,i])==F)
FCMeth<-rbind(FCMeth,cbind(BinTab[a,2],BinTab[a,i]))
}
library(lme4)
colnames(FCMeth)<-c("gene","BinaryMeth")
FCMeth[,2]<-as.numeric(FCMeth[,2])
LMRandFC<-glmer(BinaryMeth~1|gene, family=binomial(link="logit"),data=FCMeth)
GeneEstimatesFC<-ranef(LMRandFC)[[1]]
plot(density(GeneEstimatesFC[,1]), xlab="methylationEffect", main="FCFC methylation")
dev.copy(pdf, "RandomEffect_FCFC.pdf")
dev.off()



NCMeth<-data.frame()
for(i in 15:18){
a<-which(is.na(BinTab[,i])==F)
NCMeth<-rbind(NCMeth,cbind(BinTab[a,2],BinTab[a,i]))
}

NCMeth[,2]<-as.numeric(NCMeth[,2])
colnames(NCMeth)<-c("gene","BinaryMeth")
LMRandNC<-glmer(BinaryMeth~1|gene,family=binomial(link="logit"),data=NCMeth)



GeneEstimatesFC<-ranef(LMRandFC)[[1]]
GeneEstimatesNC<-ranef(LMRandNC)[[1]]


#do genes that change from FC to NC in me group change in expression?
# geneXPSN<-read.table("normalized_counts.txt", sep="\t", stringsAsFactors=F, header=T)
# 
# sampleInfo<-read.table("seqsamples.csv", sep=",",header=T)
# 
# FCFCfam<-sampleInfo[which(sampleInfo$Evol_Env=="FC_FC"),1]
# NCNCfam<-sampleInfo[which(sampleInfo$Evol_Env=="NC_NC"),1]
# 
# FCFCmean<-apply(geneXPSN[,FCFCfam],1,mean)
# NCNCmean<-apply(geneXPSN[,NCNCfam],1,mean)

MethFC<-row.names(GeneEstimatesFC)[which(GeneEstimatesFC[,1]>2.5)]
MethNC<-row.names(GeneEstimatesNC)[which(GeneEstimatesNC[,1]>2.5)]
UnMethFC<-row.names(GeneEstimatesFC)[which(GeneEstimatesFC[,1]<2)]
UnMethNC<-row.names(GeneEstimatesNC)[which(GeneEstimatesNC[,1]<2)]



#find the correlation between CpGs in methylation level within same gene and between genes
GetMePairs<-function(data_in,geneset){

data_in<-data_in[which(data_in[,2]%in%geneset==T),]

GeneCorr<-c()

for(i in 1:length(geneset)){

meanSet<-data_in[which(data_in[,2]==geneset[i]),5]


GeneCorr<-rbind(GeneCorr, expand.grid(meanSet,meanSet))


}
return(GeneCorr)
}

CorBlocks<-list()
for(i in 1:50){L<-sample(MethFC,100); CorBlocks[[i]]<-GetMePairs(MeanTab,L)}

CorReal<-c()
for(i in 1:50){CorReal<-c(CorReal, cor.test(CorBlocks[[i]][,1],CorBlocks[[i]][,2])$estimate)}

ShuffBlocks<-list()

for(i in 1:50){L<-sample(MethFC,100)
Shuff<-MeanTab
Shuff[,5]<-sample(MeanTab[,5], nrow(MeanTab))
ShuffBlocks[[i]]<-GetMePairs(Shuff,L)}
CorShuff<-c()

for(i in 1:50){CorShuff<-c(CorShuff, cor.test(ShuffBlocks[[i]][,1], ShuffBlocks[[i]][,2])$estimate)}

MethFCnoNC<-intersect(MethFC,UnMethNC)

CorRealFCnoNC<-GetMePairs(MeanTab, MethFCnoNC)

CorRealFCnoNC

colMeans(CorRealFCnoNC, na.rm=T)

######
MethNC<-intersect(MethNC,UnMethFC)

CorRealNCnoFC<-GetMePairs(MeanTab, MethNC)
######


FCNCMeth<-data.frame()
for(i in 8:11){

a<-which(is.na(BinTab[,i])==F)
FCNCMeth<-rbind(FCNCMeth,cbind(BinTab[a,2],BinTab[a,i]))
}
library(lme4)
colnames(FCNCMeth)<-c("gene","BinaryMeth")
FCNCMeth[,2]<-as.numeric(FCNCMeth[,2])
LMRandFCNC<-glmer(BinaryMeth~1|gene, family=binomial(link="logit"),data=FCNCMeth)


NCFCMeth<-data.frame()
for(i in 12:14){
a<-which(is.na(BinTab[,i])==F)
NCFCMeth<-rbind(NCFCMeth,cbind(BinTab[a,2],BinTab[a,i]))
}

NCFCMeth[,2]<-as.numeric(NCMeth[,2])
colnames(NCFCMeth)<-c("gene","BinaryMeth")
LMRandNC<-glmer(BinaryMeth~1|gene,family=binomial(link="logit"),data=NCFCMeth)



GeneEstimatesFCNC<-ranef(LMRandFCNC)[[1]]
GeneEstimatesNCFC<-ranef(LMRandNCFC)[[1]]
methFCNC<-row.names(GeneEstimatesFCNC)[which(GeneEstimatesFCNC[,1]>2.5)]
unmethFCNC<-row.names(GeneEstimatesFCNC)[which(GeneEstimatesFCNC[,1]<2)]

methFCNConly<-intersect(methFCNC,UnMethFC)
methFCFCandNC<-intersect(methFCNC,MethFC)
unmethFCNCnoFC<-intersect(unmethFCNC,MethFC)
unmethFCFCandNC<-intersect(unmethFCNC,UnMethFC)

CorRealFCNC<-GetMePairs(MeanTab,methFCNConly)
