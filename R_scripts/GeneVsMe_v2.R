# Script to identify methylated genes between conditions

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

BinTab<-CpGall

for(i in 5:ncol(CpGall)){
BinTab[,i]<-BinMeth(CpGall[,i])



}

#now we can look for CpGs that show altered methylation and test whether those within the same gene are more likely to change concordantly.

#First though we have to filter for coverage in all replicates because otherwise the tendency for high coverage over some parts of the genome will influence the output

#first comparison is between FCFC and NCNC

CoveredFN<-BinTab[which(is.na(apply(BinTab[,5:7],1,min))==F&is.na(apply(BinTab[,15:18],1,min))==F),]
UpFCFC<-which(apply(CoveredFN[,5:7],1,min)==1&apply(CoveredFN[,15:18],1,max)==0)
UpNCNC<-which(apply(CoveredFN[,5:7],1,max)==0&apply(CoveredFN[,15:18],1,min)==1)

#no sites that match this. 
#too stringent to require coverage for all sites.

#therefore alter code as follows
UpFCFC<-which(apply(BinTab[,5:7],1,mean,na.rm=T)==1&apply(BinTab[,15:18],1,mean,na.rm=T)==0)
UpNCNC<-which(apply(BinTab[,5:7],1,mean,na.rm=T)==0&apply(BinTab[,15:18],1,mean,na.rm=T)==1)

RandUpFC<-c()
RandUpNC<-c()
for(i in 1:100){
a<-sample(BinTab$gene,length(BinTab$gene),replace=F)
RandUpFC<-c(RandUpFC,length(unique(a[UpFCFC])))
RandUpNC<-c(RandUpNC,length(unique(a[UpNCNC])))
}


plot(density(RandUpFC+RandUpNC),xlim=c(1600,2200),xlab="number of genes",ylab="proportion of simulations",main="")
abline(v=2010, col="blue")
text(2100, 0.02, "real", col="blue")
text(1750, 0.02, "simulations")
text(1900, 0.02, "p<0.01\nempirical")
dev.copy(pdf, "CpG_shared_tendency_new.pdf")
dev.off()


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
geneXPSN<-read.table("normalized_counts.txt", sep="\t", stringsAsFactors=F, header=T)

sampleInfo<-read.table("seqsamples.csv", sep=",",header=T)

FCFCfam<-sampleInfo[which(sampleInfo$Evol_Env=="FC_FC"),1]
NCNCfam<-sampleInfo[which(sampleInfo$Evol_Env=="NC_NC"),1]

FCFCmean<-apply(geneXPSN[,FCFCfam],1,mean)
NCNCmean<-apply(geneXPSN[,NCNCfam],1,mean)

MethFC<-row.names(GeneEstimatesFC)[which(GeneEstimatesFC[,1]>2.5)]
MethNC<-row.names(GeneEstimatesNC)[which(GeneEstimatesNC[,1]>2.5)]
UnMethFC<-row.names(GeneEstimatesFC)[which(GeneEstimatesFC[,1]<2)]
UnMethNC<-row.names(GeneEstimatesNC)[which(GeneEstimatesNC[,1]<2)]

MethFCnoNC<-intersect(MethFC,UnMethNC)
MethNCnoFC<-intersect(MethNC,UnMethFC)
MethBoth<-intersect(MethFC,MethNC)
UnMethBoth<-intersect(UnMethFC,UnMethNC)

lookup<-read.table("gene_names.csv", sep=",", stringsAsFactors=F,fill=T,header=T)
MethFCnoNC<-cbind(MethFCnoNC,lookup[match(MethFCnoNC,lookup[,1]),2])
MethBoth<-cbind(MethBoth, lookup[match(MethBoth, lookup[,1]),2])
MethNCnoFC<-cbind(MethNCnoFC, lookup[match(MethNCnoFC,lookup[,1]),2])
UnMethBoth<-cbind(UnMethBoth, lookup[match(UnMethBoth, lookup[,1]),2])

boxplot(log2(FCFCmean[MethFCnoNC[,2]])-log2(NCNCmean[MethFCnoNC[,2]]), log2(FCFCmean[MethBoth[,2]])-log2(NCNCmean[MethBoth[,2]]),log2(FCFCmean[MethNCnoFC[,2]])-log2(NCNCmean[MethNCnoFC[,2]]),log2(FCFCmean[UnMethBoth[,2]])-log2(NCNCmean[UnMethBoth[,2]]),names=c("methFConly","methboth","methNConly","unmethBoth"),ylab="log2_meanXpsnFCFC/NCNC",notch=T, outline=F)

dev.copy(pdf, "changeInExpression_methylatedGenes.pdf")
dev.off()
#no clear change between meth both, meth FC only or meth NC only.  Interestingly though methylated genes are more highly expressed in FCFC populations than in NCNC populations, whilst the reverse is true for genes that are unmethylated in both.  

ExpVarFC<-(apply(geneXPSN[,FCFCfam],1,sd)/apply(geneXPSN[,FCFCfam],1,mean))^2
ExpVarNC<-(apply(geneXPSN[,NCNCfam],1,sd)/apply(geneXPSN[,NCNCfam],1,mean))^2

LoessFC<-loess(log2(ExpVarFC)~log2(FCFCmean))
LoessNC<-loess(log2(ExpVarNC)~log2(NCNCmean))

ResidVarFC<-residuals(LoessFC)
ResidVarNC<-residuals(LoessNC)

boxplot(ResidVarFC[MethFCnoNC[,2]]-ResidVarNC[MethFCnoNC[,2]],ResidVarFC[MethNCnoFC[,2]]-ResidVarNC[MethNCnoFC[,2]],ResidVarFC[MethBoth[,2]]-ResidVarNC[MethBoth[,2]],ResidVarFC[UnMethBoth[,2]]-ResidVarNC[UnMethBoth[,2]],notch=T,names=c("methFConly","MethNConly","methBoth","unmethBoth"), ylab="variabilityResidualFC-NC",outline=F)
wilcox.test(ResidVarFC[MethFCnoNC[,2]]-ResidVarNC[MethFCnoNC[,2]],ResidVarFC[MethNCnoFC[,2]]-ResidVarNC[MethNCnoFC[,2]],alternative="less")
wilcox.test(ResidVarFC[MethBoth[,2]]-ResidVarNC[MethBoth[,2]],ResidVarFC[UnMethBoth[,2]]-ResidVarNC[UnMethBoth[,2]],alternative="less")
text(1.5, 2.5, "p=0.03")
text(2.5, 1.5, "p=1.7e-43")
dev.copy(pdf, "variance_vsMethylation.pdf")
dev.off()


#a more sophisticated analysis is to consider genes that are EITHER environmental OR evolved separately.
#to identify environmental but reversible changes we look for genes that are up FCFC vs FCNC and up FCNC v FCFC 


EnvFCFCup<-which(apply(BinTab[,5:7],1,mean,na.rm=T)==1&apply(BinTab[,8:11],1,mean,na.rm=T)==0)
EnvFCFCdown<-which(apply(BinTab[,5:7],1,mean,na.rm=T)==0&apply(BinTab[,8:11],1,mean,na.rm=T)==1)

RandEnvUpFC<-c()
RandEnvdownNC<-c()
for(i in 1:100){
a<-sample(BinTab$gene,length(BinTab$gene),replace=F)
RandEnvUpFC<-c(RandEnvUpFC,length(unique(a[EnvFCFCup])))
RandEnvdownNC<-c(RandEnvdownNC,length(unique(a[EnvFCFCdown])))
}


plot(density(RandEnvUpFC+RandEnvdownNC),xlab="number of genes",ylab="proportion of simulations",main="",xlim=c())
abline(v=length(unique(a[EnvFCFCup]))+length(unique(a[EnvFCFCdown])), col="blue")
text(3380, 0.01, "real", col="blue")
text(3250, 0.01, "simulations", col="black")
#demonstrates some tendency for environmental genes FC to behave en bloc although not as large
dev.copy(pdf, "envFCFCvsFCNC_sharedTendency.pdf")
dev.off()

#for completeness we can now do exactly the same analysis but for NC genes that reverse

EnvNCNCup<-which(apply(BinTab[,15:18],1,mean,na.rm=T)==1&apply(BinTab[,12:14],1,mean,na.rm=T)==0)
EnvNCNCdown<-which(apply(BinTab[,15:18],1,mean,na.rm=T)==0&apply(BinTab[,12:14],1,mean,na.rm=T)==1)

RandEnvNCNCup<-c()
RandEnvNCNCdown<-c()
for(i in 1:100){
a<-sample(BinTab$gene, length(BinTab$gene),replace=F)
RandEnvNCNCup<-c(RandEnvNCNCup,length(unique(a[EnvNCNCup])))
RandEnvNCNCdown<-c(RandEnvNCNCdown,length(unique(a[EnvNCNCdown])))

}

plot(density(RandEnvNCNCup+RandEnvNCNCdown),xlab="number of genes", ylab="proportion of simulations", main="")
abline(v=length(unique(a[EnvNCNCup]))+length(unique(a[EnvNCNCdown])),col="blue")
text(1420, 0.025, "real")
text(1480, 0.025, "simulations")
#also demonstrates no tendency for environmental genes NC to behave en bloc
dev.copy(pdf, "envNCNCvsFCNC_sharedTendency.pdf")
dev.off()

#can say then that evolution affects methylation of some genes but environment has little influence. 

#could argue then that there is no coherent justification for identifying genes that show differential methylation env NC or env FC but could do this anyway.  


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
FCNCfam<-sampleInfo[which(sampleInfo$Evol_Env=="FC_NC"),1]
FCNCmean<-apply(geneXPSN[,FCNCfam],1,mean)

methFCNConly<-cbind(methFCNConly, lookup[match(methFCNConly,lookup[,1]),2])
methFCFCandNC<-cbind(methFCFCandNC, lookup[match(methFCFCandNC, lookup[,1]),2])
unmethFCNCnoFC<-cbind(unmethFCNCnoFC,lookup[match(unmethFCNCnoFC,lookup[,1]),2])
unmethFCFCandNC<-cbind(unmethFCFCandNC,lookup[match(unmethFCFCandNC,lookup[,1]),2])


boxplot(log2(FCNCmean[methFCNConly[,2]])-log2(FCFCmean[methFCNConly[,2]]),log2(FCNCmean[methFCFCandNC[,2]])-log2(FCFCmean[methFCFCandNC[,2]]),log2(FCNCmean[unmethFCNCnoFC[,2]])-log2(FCFCmean[unmethFCNCnoFC[,2]]),log2(FCNCmean[unmethFCFCandNC[,2]])-log2(FCFCmean[unmethFCFCandNC[,2]]),names=c("methFCNConly","methboth","unmethFCNConly","unmethBoth"),ylab="log2(FCNC/FCFC)",outline=F,notch=T)
#this is an odd result.  Genes that are methylated either in FCFC or FCNC tend to have lower expression in FCNC, not true for genes that are unmethylated in both. This suggests that the tendency for methylated genes to be higher in expression in FCFC is environmental.  Could fit with the idea of hk genes being more highly expressed in FCFC, which is to be investigated later.

dev.copy(pdf, "methylation_effectFCNCvFCFC_expression.pdf")
dev.off()
#what about variance? 
#to look at this we calculate the residuals of the loess fit cv2 against mean, as above.
ExpVarFCNC<-(apply(geneXPSN[,FCNCfam],1,sd)/apply(geneXPSN[,FCNCfam],1,mean))^2
LoessFCNC<-loess(log2(ExpVarFCNC)~log2(FCNCmean))
ResidVarFCNC<-residuals(LoessFCNC)
boxplot(ResidVarFCNC[methFCNConly[,2]]-ResidVarFC[methFCNConly[,2]],ResidVarFCNC[unmethFCNCnoFC[,2]]-ResidVarFC[unmethFCNCnoFC[,2]],ResidVarFCNC[methFCFCandNC[,2]]-ResidVarFC[methFCFCandNC[,2]],ResidVarFCNC[unmethFCFCandNC[,2]]-ResidVarFC[unmethFCFCandNC[,2]],outline=F,notch=T, names=c("methFCNConly","methFCFConly","methBoth","unmethBoth"),ylab="variability FCNC-FCFC")
#this again fits with the idea that methylation tends to reduce variability and loss of methylation increase it. Interestingly there is more variability in expression in the FCNC relative to the FCFC but specifically for genes that are methylated in either. 
dev.copy(pdf, "methylation_effectonvariability_envFCFCvFCNC.pdf")
dev.off()


#another fundamental question is whether methylated genes are more variable or not
MethFC<-cbind(MethFC, lookup[match(MethFC,lookup[,1]),2])
UnMethFC<-cbind(UnMethFC,lookup[match(UnMethFC,lookup[,1]),2])
MethNC<-cbind(MethNC,lookup[match(MethNC, lookup[,1]),2])

UnMethNC<-cbind(UnMethNC, lookup[match(UnMethNC,lookup[,1]),2])
boxplot(ResidVarFC[MethFC[,2]],ResidVarFC[UnMethFC[,2]],ResidVarNC[MethNC[,2]],ResidVarNC[UnMethNC[,2]],names=c("methFC","unmethFC","methNC","unmethNC"),notch=T, outline=F, ylab="residuals_CV2vmean")
dev.copy(pdf, "effect_of_methylation_on_variability_allgenes.pdf")
dev.off()
#can also be plotted as cv vs mean graph

plot(log2(FCFCmean),log2(ExpVarFC),pch=16,ylab="expression variability (CV2)",xlab="mean")
points(log2(FCFCmean[MethFC[,2]]),log2(ExpVarFC[MethFC[,2]]), pch=16,col="red")
legend("bottomleft", c("unmethylated","methylated"), col=c("black","red"), pch=16)
dev.copy(pdf, "effect_of_methylation_on_Cv2_scatter.pdf")
dev.off()

#save.image("data_analysis_0910.Rdata")
snp<-read.table("snps/genes_intersect_snp.bed", sep="\t", stringsAsFactors=F)

colNames<-read.table("snps/column_names.csv", sep=",", stringsAsFactors=F,header=T)

colnames(snp)[13:54]<-paste(colNames[,3],colNames[,4], sep=":")

FCFCsnp<-grep("freq:FCFC", colnames(snp))
NCNCsnp<-grep("freq:NCNC", colnames(snp))

MeanSnpFreqFC<-c()
for(i in 1:length(FCFCsnp)){MeanSnpFreqFC<-cbind(MeanSnpFreqFC, as.numeric(unlist(snp[,FCFCsnp[i]])))}
MeanSnpFreqNC<-c()
for(i in 1:length(NCNCsnp)){MeanSnpFreqNC<-cbind(MeanSnpFreqNC, as.numeric(unlist(snp[,NCNCsnp[i]])))}
MeanSnpFreqFC<-apply(MeanSnpFreqFC,1,mean,na.rm=T)
MeanSnpFreqNC<-apply(MeanSnpFreqNC,1,mean,na.rm=T)


FCvsNCsnp<-which(MeanSnpFreqFC>0.8&MeanSnpFreqNC<=0.2)
NCvsFCsnp<-which(MeanSnpFreqFC<0.2&MeanSnpFreqNC>=0.8)

GenesSnpFCvsNC<-read.table(text=snp[FCvsNCsnp,7], sep=";")[,1]
GenesSnpNCvsFC<-read.table(text=snp[NCvsFCsnp,7],sep=";")[,2]

GenesMethFC_withSnp<-unique(c(intersect(MethFCnoNC[,1],GenesSnpFCvsNC),intersect(MethFCnoNC[,1],GenesSnpNCvsFC)))
GenesMethNC_withSnp<-unique(c(intersect(MethNCnoFC[,1],GenesSnpFCvsNC),intersect(MethNCnoFC[,1],GenesSnpNCvsFC)))
GenesMethBoth_withSnp<-unique(c(intersect(MethBoth[,1],GenesSnpFCvsNC),intersect(MethBoth[,1],GenesSnpNCvsFC)))
GenesUnmethBoth_withSnp<-unique(c(intersect(UnMethBoth[,1],GenesSnpFCvsNC)),intersect(UnMethBoth[,1],GenesSnpNCvsFC))

Snp_table<-c(length(GenesMethFC_withSnp),length(GenesMethNC_withSnp),length(GenesMethBoth_withSnp),length(GenesUnmethBoth_withSnp))
Snp_table<-rbind(Snp_table, c(nrow(MethFCnoNC),nrow(MethNCnoFC),nrow(MethBoth),nrow(UnMethBoth)))

barplot(Snp_table[1,]/Snp_table[2,]*100,ylim=c(0,4),names=c("FCFCmeth_only","NCNCmeth_only","Bothmeth","bothUnmeth"),las=2, ylab="%genes with differential snps")
fisher.test(rbind(Snp_table[1,3:4],Snp_table[2,3:4]-Snp_table[1,3:4]))

text(3.5,3.9, "Padj=2e-15\nmeth v unmeth")


#table of gene estimates 
AllGenes<-c(row.names(GeneEstimatesFC),row.names(GeneEstimatesNC),row.names(GeneEstimatesNCFC),row.names(GeneEstimatesFCNC))
AllGenes<-unique(AllGenes)
GeneTab_all<-matrix(0, ncol=4, nrow=length(AllGenes))
row.names(GeneTab_all)<-AllGenes
GeneEstimateList<-list(GeneEstimatesFC,GeneEstimatesNC,GeneEstimatesNCFC,GeneEstimatesFCNC)
names(GeneEstimateList)<-c("FCFC","NCNC","NCFC","FCNC")
for(i in 1:4){
q<-match(row.names(GeneTab_all),row.names(GeneEstimateList[[i]]))
GeneTab_all[,i]<-GeneEstimateList[[i]][q,1]



}
colnames(GeneTab_all)<-names(GeneEstimateList)
GeneTab_allBin<-GeneTab_all
for(i in 1:4){GeneTab_allBin[which(GeneTab_all[,i]>2.5),i]<-1
GeneTab_allBin[which(GeneTab_all[,i]<2),i]<-0
GeneTab_allBin[-which(GeneTab_allBin[,i]==0|GeneTab_allBin[,i]==1),i]<-(-1)
}

