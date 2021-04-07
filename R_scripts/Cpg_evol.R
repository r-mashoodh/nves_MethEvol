## methylkit and cpg methylation

library(methylKit)
library(tidyverse)

# Create a list of files of Me data
file.list <- list("../Evol_Methylation/FCFC_01_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/FCFC_02_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/FCFC_03_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/FCNC_01_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/FCNC_02_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/FCNC_03_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/FCNC_04_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/NCFC_01_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/NCFC_02_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/NCFC_03_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/NCNC_01_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/NCNC_02_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/NCNC_03_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "../Evol_Methylation/NCNC_04_1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz.bismark.zero.cov")

# code evolved condition as group (FC = 0, NC = 1)
# add covariate?
raw_data <- methRead(file.list,
                     sample.id = list("FCFC01", "FCFC02", "FCFC03",
                                      "FCNC01", "FCNC02", "FCNC03", "FCNC04",
                                      "NCFC01", "NCFC02", "NCFC03",
                                      "NCNC01", "NCNC02", "NCNC03", "NCNC04"),
                     treatment = c(0, 0, 0,
                                   0, 0, 0, 0,
                                   1, 1, 1,
                                   1, 1, 1, 1),
                     assembly="Nves",
                     context="CpG",
                     pipeline = "bismarkCoverage")

# get some descriptive stats on samples

# # gives you summary and percentiles
# getMethylationStats(raw_data[[2]],plot=FALSE,both.strands=FALSE)
#
# # plot of percent methylation
# # most is less than 20% methylation
# getMethylationStats(raw_data[[2]],plot=TRUE,both.strands=FALSE)
#
# # coverage plot
# getCoverageStats(raw_data[[2]],plot=TRUE,both.strands=FALSE)


# Filter data for outliers and coverage and normalise counts
filtered_data <- filterByCoverage(raw_data,
                                  lo.count=10,
                                  lo.perc=NULL,
                                  hi.count=NULL,
                                  hi.perc=99.9)


meth_all_data <- methylKit::unite(filtered_data, min.per.group=5L)
# 
# 
df_meth_all <- getData(meth_all_data)

# Ugly way of changing every NA coverage value to 100
df_meth_all[,c(5,8,11,14,17,20,23,26,29,32,35,38,41,44)][is.na(df_meth_all[,c(5,8,11,14,17,20,23,26,29,32,35,38,41,44)])] <- 100
# Ugly way of changing every numCs value to 1
df_meth_all[,c(6,9,12,15,18,21,24,27,30,33,36,39,42,45)][is.na(df_meth_all[,c(6,9,12,15,18,21,24,27,30,33,36,39,42,45)])] <- 1
# Using these values means the NA site is not called as methylated
# This is a lazy way of doing it and not sensible if you have 100s of samples

a <- df_meth_all[,5:6]
b <- df_meth_all[,8:9]

c <- df_meth_all[,11:12]
d <- df_meth_all[,14:15]

e <- df_meth_all[,17:18]
f <- df_meth_all[,20:21]

g <- df_meth_all[,23:24]
h <- df_meth_all[,26:27]

i <- df_meth_all[,29:30]
j <- df_meth_all[,32:33]

k <- df_meth_all[35:36]
l <- df_meth_all[38:39]

m <- df_meth_all[41:42]
n <- df_meth_all[44:45]

# NOTE: p shouold be the average non-conversion rate (proportion of methylated Cs compared to non-meth Cs)
# So if 1000 methylated Cs compared to 200,000 T's then 1000/200,000 = 0.005
# for a paper: 'the success probability is the non-conversion rate'
bt <- function(a, b, p = 0.005) {binom.test(a, b, 0.005, alternative="greater")$p.value}
allrows <- data.frame(CT=numeric(),Ccount=numeric(),pVal=numeric(),FDR=numeric(),row=numeric())

for (df in list(a,b,c,d,e,f,g,h,j,k,l,m,n)) {
  colnames(df) <- c("CT", "Ccount")
  df$pVal <- mapply(bt, df$Ccount, df$CT)
  df$FDR <- p.adjust(df$pVal, method = "BH", n = length(df$pVal))
  df$row <- as.numeric(rownames(df))
  dfmeth <- subset(df, df$FDR < 0.1)
  allrows <- rbind(allrows, dfmeth)
}

meth_positions <- as.vector(as.numeric(unique(allrows$row))) 
length(meth_positions) # 7400

## create a methyl object with just the 'methylated' CpGs
subset_methBase <- methylKit::select(meth_all_data, meth_positions)

#############################################################################################
## get methylated CpGs and look at distribution
MethylatedCpGs <- methylKit::getData(subset_methBase)

MethylatedCpGs %>% 
  select(chr, start, end) %>% 
  write.table(., "MethylatedCpGs.bed", sep="\t", quote=F, row.names = F, col.names = F)

# in terminal: bedtools intersect -a MethylatedCpGs.bed -b Nves.filtered.gff -wao > MethylatedCpGs_intersect.bed

m <- read.table("MethylatedCpGs_intersect_filtered.bed", sep="\t", na.strings = ".") %>% 
  rename(chr = V1,
         start = V2,
         end = V3)

m %>% 
  ggplot(., aes(x=reorder(V6, V6, function(x)-length(x)), fill=V6)) +
  geom_bar(alpha=0.6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x=element_blank())

m %>% 
  group_by(V6) %>% 
  summarise(n=n(), 
            prop=n()/7400)

# meth <- data.frame(methylKit::percMethylation(subset_methBase)) %>% 
#   mutate(mean_CpG = rowMeans(., na.rm=T)) %>% 
#   select(mean_CpG)
# 
# MethylatedCpGs <- cbind(MethylatedCpGs, meth)
# 
# m %>% 
#   left_join(MethylatedCpGs) %>% 
#   ggplot(., aes(x=reorder(V6, V6, function(x)-length(x)), y=log2(mean_CpG), color=V6)) +
#   geom_violin() +
#   stat_summary(fun.y=mean, geom="point", shape=23, size=2)

############################################################################################

## Differential methylation for 'Evol'

meth_evol <- methylKit::calculateDiffMeth(subset_methBase, slim=FALSE)
meth_evol <- methylKit::getMethylDiff(meth_evol, qvalue=0.05, difference = 0)
meth_evol
### create a bedfile of coordinates to intersect with gff
methylKit::getData(meth_evol) %>% 
  dplyr::select(chr, start, end) %>%
  write.table(., "meth_evol_sites.bed", quote=F, sep="\t", col.names = F, row.names = F)

system("bedtools intersect -a meth_evol_sites.bed -b Nves.filtered.gff -wao > meth_evol_intersect.bed")

evol_intersect <- read.table("meth_evol_intersect.bed", sep="\t", na.strings = ".")

evol_intersect %>% 
  ggplot(., aes(x=reorder(V6, V6, function(x)-length(x)))) +
  geom_bar(alpha=0.6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x=element_blank())

evol_intersect %>% 
  group_by(V6) %>% 
  summarise(n=n(), 
            prop=n()/536)

######################################################################################################
## Reorganise the methylkit object to test for effects of 'env' 

meth.reorg <-reorganize(subset_methBase, sample.ids = c("FCFC01", "FCFC02", "FCFC03",
                                                        "FCNC01", "FCNC02", "FCNC03", "FCNC04",
                                                        "NCFC01", "NCFC02", "NCFC03",
                                                        "NCNC01", "NCNC02", "NCNC03", "NCNC04"),
                        treatment = c(0, 0, 0,
                                      1, 1, 1, 1,
                                      0, 0, 0,
                                      1, 1, 1, 1))


meth_env <-  methylKit::calculateDiffMeth(meth.reorg, slim=F)
meth_env <- methylKit::getMethylDiff(meth_env, qvalue=0.05, difference = 0)
meth_env ## 55 sites diff methylated

### create a bedfile of coordinates to intersect with gff
methylKit::getData(meth_env) %>% 
  dplyr::select(chr, start, end) %>%
  write.table(., "meth_env_sites.bed", quote=F, sep="\t", col.names = F, row.names = F)

system("bedtools intersect -a meth_env_sites.bed -b Nves.filtered.gff -wao > meth_env_intersect.bed")

env_intersect <- read.table("meth_env_intersect.bed", sep="\t", na.strings = ".")

env_intersect %>% 
  ggplot(., aes(x=reorder(V6, V6, function(x)-length(x)))) +
  geom_bar(alpha=0.6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x=element_blank())

env_intersect %>% 
  group_by(V6) %>% 
  summarise(n=n(), 
            prop=n()/55)


#################################################################
# quick plot of proportions

prop <- evol_intersect %>% 
  group_by(V6) %>% 
  summarise(n=n(), 
            prop=n()/536) %>% 
  bind_rows(env_intersect %>% 
              group_by(V6) %>% 
              summarise(n=n(), 
                        prop=n()/55)) %>% 
  rename(feature = V6) %>% 
  mutate(exp = c(rep("evol", 5), rep("env", 5)),
         feature = factor(.$feature, levels = c("gene", "CDS", "intron", "five_prime_UTR", "TE", "NA")))

ggplot(prop, aes(x=feature, y=prop, fill=exp)) +
  geom_bar(stat="identity", position="dodge") +
  ggtitle("all p's")

########################################################################################################
## Create a dataframe to do Fishers tests for enrichment

## what mCpG are in genes 
m_genes <- m %>% 
  filter(V6 == "gene") %>% 
  mutate(siteID = paste(V1, V2, sep="_"))

m_TEs <- m %>% 
  filter(V6 == "TE") %>% 
  mutate(siteID = paste(V1, V2, sep="_"))

m_CDS <- m %>% 
  filter(V6 == "CDS") %>% 
  mutate(siteID = paste(V1, V2, sep="_"))

m_intron <- m %>% 
  filter(V6 == "intron") %>% 
  mutate(siteID = paste(V1, V2, sep="_"))

m_5UTR <- m %>% 
  filter(V6 == "five_prime_UTR") %>% 
  mutate(siteID = paste(V1, V2, sep="_"))

# m_upstream <- m %>% 
#   filter(V6 == "gene_upstream_flank") %>% 
#   mutate(siteID = paste(V1, V2, sep="_"))
# 
# m_downstream <- m %>% 
#   filter(V6 == "gene_downstream_flank") %>% 
#   mutate(siteID = paste(V1, V2, sep="_"))

m_NA <- m %>% 
  filter(is.na(V6)) %>% 
  mutate(siteID = paste(V1, V2, sep="_"))

## what CpGs are DM?
dm_evol <- methylKit::getData(meth_evol) %>% 
  mutate(siteID = paste(chr, start, sep="_"))

dm_env <- methylKit::getData(meth_env) %>% 
  mutate(siteID = paste(chr, start, sep="_"))


# get full list of methylated Cpgs
# mark the genes
Intersect.df <- MethylatedCpGs %>%  
  dplyr::select(chr, start, end) %>% 
  mutate(siteID = paste(chr, start, sep="_")) %>% 
  mutate(dm_evol = ifelse(siteID %in% dm_evol$siteID, 1, 0),
         dm_env = ifelse(siteID %in% dm_env$siteID, 1, 0),
         mGene= ifelse(siteID %in% m_genes$siteID, 1, 0),
         mTE = ifelse(siteID %in% m_TEs$siteID, 1, 0),
         mCDS = ifelse(siteID %in% m_CDS$siteID, 1, 0),
         mIntron = ifelse(siteID %in% m_intron$siteID, 1, 0),
         m5UTR = ifelse(siteID %in% m_5UTR$siteID, 1, 0),
         #mUp = ifelse(siteID %in% m_upstream$siteID, 1,0),
         #mDown = ifelse(siteID %in% m_upstream$siteID, 1, 0),
         mNA = ifelse(siteID %in% m_NA$siteID, 1, 0))

#### get Fishers test results 
Fishers.evol <- apply(Intersect.df[,7:12], 2, fisher.test, y = Intersect.df$dm_evol)
Fishers.env <- apply(Intersect.df[,7:12], 2, fisher.test, y = Intersect.df$dm_env)

#############################
## plotting

evol.tmp <- as.data.frame(unlist(lapply(Fishers.evol, function(x) x$conf.int))) 
colnames(evol.tmp)[1] <- c("conf")

evol.results <- evol.tmp %>% 
  rownames_to_column(var = "rowname") %>% 
  mutate(int = substrRight(.$rowname, 1)) %>% 
  mutate(int = ifelse(int == 1, "lwr", "upr")) %>% 
  mutate(feature = gsub("[[:digit:]]","",.$rowname),
         value = as.numeric(.$conf)) %>% 
  select(-rowname, -conf) %>% 
  pivot_wider(names_from=int) %>% 
  mutate(odds = unlist(lapply(Fishers.evol, function(x) x$estimate)),
         pval = unlist(lapply(Fishers.evol, function(x) x$p.value))) %>% 
  mutate(padj = p.adjust(.$pval, method="BH", n=6))

########################################################

env.tmp <- as.data.frame(unlist(lapply(Fishers.env, function(x) x$conf.int))) 
colnames(env.tmp)[1] <- c("conf")

env.results <- env.tmp %>% 
  rownames_to_column(var = "rowname") %>% 
  mutate(int = substrRight(.$rowname, 1)) %>% 
  mutate(int = ifelse(int == 1, "lwr", "upr")) %>% 
  mutate(feature = gsub("[[:digit:]]","",.$rowname),
         value = as.numeric(.$conf)) %>% 
  select(-rowname, -conf) %>% 
  pivot_wider(names_from=int) %>% 
  mutate(odds = unlist(lapply(Fishers.env, function(x) x$estimate)),
         pval = unlist(lapply(Fishers.env, function(x) x$p.value))) %>% 
  mutate(padj = p.adjust(.$pval, method="BH", n=6))


full <- evol.results %>% 
  bind_rows(env.results) %>% 
  mutate(exp = c(rep("evol",6), rep("env",6)))

full$feature <- factor(full$feature, levels = c("mGene", "mCDS", "mIntron", "mUTR", "mTE", "mNA"))


ggplot(full, aes(x=feature, y=odds, color=feature)) +
  geom_point() +
  geom_errorbar(aes(ymin=lwr,ymax=upr,color=feature,width=0.2)) +
  geom_hline(yintercept = 1, color="grey") +
  facet_wrap(.~exp) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  ylim(0,6) +
  ggtitle("all sig p's")


# library(gridExtra)
# 
# grid.arrange(ggplot(full100, aes(x=feature, y=odds, fill=exp)) +
#                geom_bar(stat="identity", position="dodge") +
#                ylim(0,3) +
#                ggtitle("top 100"),
#              ggplot(full, aes(x=feature, y=odds, fill=exp)) +
#                geom_bar(stat="identity", position="dodge") +
#                ylim(0,3) +
#                ggtitle("all sig"))


### combine and then do fishers test

# dm_evol and dm_env ... to combine just add?

Intersect.all <- Intersect.df %>% 
  mutate(dm_all = ifelse(dm_evol == 1 | dm_env==1, 1, 0))
  

Fishers.all <- apply(Intersect.all[,7:12], 2, fisher.test, y = Intersect.all$dm_all)
Fishers.all

all.tmp <- as.data.frame(unlist(lapply(Fishers.all, function(x) x$conf.int))) 
colnames(all.tmp)[1] <- c("conf")

all.results <- all.tmp %>% 
  rownames_to_column(var = "rowname") %>% 
  mutate(int = substrRight(.$rowname, 1)) %>% 
  mutate(int = ifelse(int == 1, "lwr", "upr")) %>% 
  mutate(feature = gsub("[[:digit:]]","",.$rowname),
         value = as.numeric(.$conf)) %>% 
  select(-rowname, -conf) %>% 
  pivot_wider(names_from=int) %>% 
  mutate(odds = unlist(lapply(Fishers.all, function(x) x$estimate)),
         pval = unlist(lapply(Fishers.all, function(x) x$p.value))) %>% 
  mutate(padj = p.adjust(.$pval, method="BH", n=6))

all.results$feature <- factor(all.results$feature, levels = c("mGene", "mCDS", "mIntron", "mTE", "mUTR", "mNA"))

###### PLOTS for PAPER #####

colors <- c("#440154FF", "#3b528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF", "grey")


cpg_dist <- m %>% 
  group_by(V6) %>% 
  summarise(n=n(), 
            prop=(n()/7400)*100) %>% 
  mutate(V6 = factor(V6, levels = c("gene", "CDS", "intron", "TE", "five_prime_UTR", "NA"))) %>% 
  ggplot(., aes(x = V6, y = prop, fill = V6)) + 
  geom_bar(stat = "identity", alpha=0.7) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x=element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = colors, 
                    na.value = "grey")

sig_cpg <- ggplot(all.results, aes(x=feature, y=odds, color=feature)) +
  geom_point() +
  geom_errorbar(aes(ymin=lwr,ymax=upr,color=feature,width=0.2)) +
  geom_hline(yintercept = 1, color="grey") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none") +
  scale_color_manual(values = colors)

(cpg_dist / sig_cpg)

density.table <- m %>% 
  group_by(V6) %>% 
  summarise(n=n(), 
            prop=(n()/7400)*100)

chisq <- chisq.test(density.table$n)
chisq


library(patchwork)

