# Feb 6 2020

# Load required libraries
library(tidyverse)
library(DESeq2)
library(knitr)
library(EnhancedVolcano)
library(eulerr)
#library(apeglm)
#library(RColorBrewer)
#library(Glimma)
#library(gplots)


# Read in counts table
genes.matrix <- read.csv("gene_count_matrix.csv")

# geneID to rownames
rownames(genes.matrix) <- genes.matrix$gene_id
genes.matrix <- genes.matrix %>% select(-gene_id)


## Create sample info file 
## get groups extracted
seqsamples <- read.csv("seqsamples.csv")
sample_info <- seqsamples %>% 
  dplyr::select(Family,Pop_Rep,Evolved,Group,Env,Evol_Env)

# renaming some groups
sample_info$rep1 <- substr(sample_info$Pop_Rep, 2, 2)
sample_info$rep.mean.center = ifelse(sample_info$rep1 == 1, -0.5, 0.5) 
sample_info$rep = ifelse(sample_info$rep1 == 1, "B1", "B2")

# check matched order
seIdx <- match(colnames(matrix), sample_info$Family) ### is this working?

# check colnames/rownames, they need to match samples
# should be true
all(colnames(genes.matrix) == sample_info$Family)

genes.matrix <- round(genes.matrix,0)

dds <- DESeqDataSetFromMatrix(countData = genes.matrix,
                              colData = sample_info,
                              design = ~ Evol_Env + rep.mean.center)

normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="normalized_counts.txt", sep="\t", quote=F)


# pre-filtering
keep <- rowSums(counts(dds)) >= 15
dds <- dds[keep,]

# estimate size factors, for normalisation
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

# check result names
resultsNames(dds)

# Check mean center coefficient
block.res = results(dds, name=c("rep.mean.center"))
block.res.shrink = lfcShrink(dds, res=block.res,
                             coef=c("rep.mean.center"),
                             type="ashr")

# effect of block
# 34 genes show block effect
block.res <- block.res.shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  arrange(desc(log2FoldChange))




# 1. Effect of NC-Env in FC pop
# FCFC vs FCNC

FCenv = results(dds, contrast=c("Evol_Env","FC_NC","FC_FC"))
FCenv.shrink = lfcShrink(dds, res=FCenv,
                         contrast=c("Evol_Env","FC_NC","FC_FC"),
                         type="ashr")

############### 
# list of all genes for randomization

all.genes.list <- FCenv.shrink %>%  
  data.frame() %>% 
  rownames_to_column(var="gene") %>% 
  select(gene)

write.table(all.genes.list, "all_genes_list.txt", col.names = F, 
            row.names = F, quote = F)
#####

## 2 fold change
FCenv.2fold <- FCenv.shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>% 
  arrange(desc(log2FoldChange))

write.csv(FCenv.2fold, "FCenv_2fold_deseq.csv", quote=F)

#up/down count
FCenv.2fold %>% 
  summarise(down = sum(log2FoldChange < 0), up = sum(log2FoldChange > 0))

## write all genes and their p-values to a file

all_FCenv <- FCenv.shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene")

write.csv(all_FCenv, "all_FCenv_deseq2_results.csv", quote=F)

## 3 fold change
# FCenv.3fold <- FCenv.shrink %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble() %>% 
#   filter(padj < 0.05 & abs(log2FoldChange) >= 1.585) %>% 
#   arrange(desc(log2FoldChange))
# 
FCenv.nofold <- FCenv.shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  filter(padj < 0.05) %>%
  arrange(desc(log2FoldChange))
# 
# #up/down count
# FCenv.3fold %>% 
#   summarise(down = sum(log2FoldChange < 0), up = sum(log2FoldChange > 0))


EnhancedVolcano(FCenv.shrink,
                lab = rownames(FCenv.shrink),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1)

## 2. Evol: Effect of 30 Gen of NC (NCEnv + unique Evol)
# FCFC vs NCNC

Evol = results(dds, contrast=c("Evol_Env","NC_NC","FC_FC"))
Evol.shrink = lfcShrink(dds, res=Evol,
                        contrast=c("Evol_Env","NC_NC","FC_FC"),
                        type="ashr")

Evol.2fold <- Evol.shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>% 
  arrange(desc(log2FoldChange))

write.csv(Evol.2fold, "Evol_2fold_deseq2.csv", quote=F)


all_Evol <- Evol.shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene")

write.csv(all_Evol, "all_Evol_deseq2_results.csv", quote=F)

# up down summary of Evol
Evol.2fold %>% 
  summarise(down = sum(log2FoldChange < 0), up = sum(log2FoldChange > 0))

EnhancedVolcano(Evol.shrink,
                lab = rownames(FCenv.shrink),
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1)

Evol.nofold <- Evol.shrink %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  filter(padj < 0.05) %>%
  arrange(desc(log2FoldChange))




# common genes with FCenv
nrow(FCenv.2fold)
length(intersect(FCenv.2fold$gene, Evol.2fold$gene))
length(setdiff(FCenv.2fold$gene, Evol.2fold$gene))
length(setdiff(Evol.2fold$gene, FCenv.2fold$gene)) #82 unique evolving genes   

### 3-fold ###
# Evol.3fold <- Evol.shrink %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble() %>% 
#   filter(padj < 0.05 & abs(log2FoldChange) >= 1.585) %>% 
#   arrange(desc(log2FoldChange))
# 
# # up down summary of Evol
# Evol.3fold %>% 
#   summarise(down = sum(log2FoldChange < 0), up = sum(log2FoldChange > 0))




## 3. NCenv diffs across pops
# 
# pops.NCenv = results(dds, contrast=c("Evol_Env","NC_NC","FC_NC"))
# pops.NCenv = lfcShrink(dds, res=pops.NCenv,
#                        contrast=c("Evol_Env","NC_NC","FC_NC"),
#                        type="ashr")
# 
# pops.NCenv <- pops.NCenv %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble() %>% 
#   filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>% 
#   arrange(desc(log2FoldChange))
# 
# 
# ## 4. FCenv diffs across pops
# pops.FCenv = results(dds, contrast=c("Evol_Env","NC_FC","FC_FC"))
# pops.FCenv = lfcShrink(dds, res=pops.FCenv,
#                        contrast=c("Evol_Env","NC_FC","FC_FC"),
#                        type="ashr")
# 
# pops.FCenv <- pops.FCenv %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble() %>% 
#   filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>% 
#   arrange(desc(log2FoldChange))


## 5. Plasticity back to FC env in NC

NCpop.FCenv = results(dds, contrast=c("Evol_Env","NC_NC","NC_FC"))
NCpop.FCenv = lfcShrink(dds, res=NCpop.FCenv,
                        contrast=c("Evol_Env","NC_NC","NC_FC"),
                        type="ashr")

NCpop.FCenv.2fold <- NCpop.FCenv  %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05 & abs(log2FoldChange) >= 1) %>% 
  arrange(desc(log2FoldChange))

NCpop.nofold <- NCpop.FCenv %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  as_tibble() %>%
  filter(padj < 0.05) %>%
  arrange(desc(log2FoldChange))

## write all genes and their p-values to a file

NC_plasticity <- NCpop.FCenv  %>%
  data.frame() %>%
  rownames_to_column(var="gene")

write.csv(NC_plasticity, "all_NCenv_deseq2_results.csv", quote=F)

## write some gene lists
# write gene lists 
# write.table(FCenv.2fold$gene, "FCenv2fold.txt", quote=FALSE, row.names = F, col.names = F)
# write.table(Evol.2fold$gene, "Evol2fold.txt", quote=FALSE, row.names = F, col.names = F)
# write.table(FCenv.3fold$gene, "FCenv3fold.txt", quote=FALSE, row.names = F, col.names = F)
# write.table(Evol.3fold$gene, "Evol3fold.txt", quote=FALSE, row.names = F, col.names = F)

# FCenv.3fold.unique <- FCenv.3fold[!FCenv.3fold$gene %in% Evol.3fold$gene,]
# Evol.3fold.unique <- Evol.3fold[!Evol.3fold$gene %in% FCenv.3fold$gene,]


#### overlap ####
## Work out how much overlap between FCenv (NC-induced plasticity) vs Evolving
nrow(FCenv.2fold)
nrow(Evol.2fold)
length(intersect(FCenv.2fold$gene, Evol.2fold$gene))
length(setdiff(FCenv.2fold$gene, Evol.2fold$gene))
length(setdiff(Evol.2fold$gene, FCenv.2fold$gene))

# ## Work out how much overlap between FCenv (NC-induced plasticity) vs Evolving
# nrow(FCenv.3fold)
# nrow(Evol.3fold)
# length(intersect(FCenv.3fold$gene, Evol.3fold$gene))
# length(setdiff(FCenv.3fold$gene, Evol.3fold$gene))
# length(setdiff(Evol.3fold$gene, FCenv.3fold$gene))

fit <- euler(c("A" = length(FCenv.2fold), "B" = length(Evol.2fold), 
               "A&B" = length(intersect(FCenv.2fold$gene, Evol.2fold$gene))), 
             shape = "circle")
             


fit$stress
fit$diagError
plot(fit)


venn.diagram(list(B = 1:1800, A = 1571:2020), fill = c("lightblue", "green"), 
             alpha = c(0.5, 0.5), lwd =0, "venn_diagram.tiff")



### Comparing reaction norms 

# Of the genes that are overlapping, are they concordantly expressed
FC_stat <- data.frame(FCenv.nofold$gene, FCenv.nofold$log2FoldChange)
colnames(FC_stat) <- c("gene", "FC.foldchange")
NC_stat <- data.frame(NCpop.nofold$gene, NCpop.nofold$log2FoldChange)
colnames(NC_stat) <- c("gene", "NC.foldchange")
log.compare <- inner_join(NC_stat, FC_stat)

sign(stat.compare.evolving$FC.foldchange) == sign(stat.compare.evolving$Evol.foldchange)

ggplot(log.compare, aes(x=FC.foldchange, y=NC.foldchange)) +
  geom_point(size=2, alpha=0.5) +
  geom_smooth(method=lm, se=TRUE) +
  geom_abline(intercept = 0, slope = 1, color="red") +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed")

library(lme4)
library(lmerTest)

env.lm <- lm(NC.foldchange ~ FC.foldchange, data=log.compare)
summary(env.lm)

r1 <- lm(NC.foldchange ~ 1 + offset(FC.foldchange), data=log.compare)
anova(env.lm, r1)

DEG.log.compare <- log.compare %>% 
  filter(gene %in% Evol.2fold$gene)

#### bar plots of updown regulation
fc.up <- FCenv.2fold %>% 
  group_by(down = log2FoldChange < 0, up = log2FoldChange > 1) %>% 
  summarise(n=n()) %>% 
  pull(n)

evol.up <- Evol.2fold %>% 
  group_by(down = log2FoldChange < 0, up = log2FoldChange > 1) %>% 
  summarise(n=n()) %>% 
  pull(n)

pops.up <- pops.FCenv %>% 
  group_by(down = log2FoldChange < 0, up = log2FoldChange > 1) %>% 
  summarise(n=n()) %>% 
  pull(n)

up.down <- data.frame(rbind(fc.up, evol.up, pops.up))
colnames(up.down) <- c("up", "down")
up.down$contrast <- row.names(up.down)

DE_long <- gather(up.down, key=DE, value=sum, 1:2)
DE_long$contrast<- factor(DE_long$contrast, levels =  c("evol.up", "fc.up", "pops.up"))

ggplot(data=DE_long, aes(x=contrast, y=sum, fill=DE)) +
  geom_bar(stat="identity", width = 0.5) +
  ylim(0, 700) +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values=c("#9ecae1", "#2171b5"))

# ggplot(data=DE_long, aes(x=contrast, y=sum, fill=DE)) +
#   geom_bar(stat="identity", width = 0.5) +
#   ylim(0, 700) +
#   theme_classic(base_size = 26) +
#   coord_flip() +
#   scale_fill_manual(values=c("#9ecae1", "#2171b5")) +
#   theme(axis.text.x = element_text(size=20))

## venn diagram

#dddddddd


