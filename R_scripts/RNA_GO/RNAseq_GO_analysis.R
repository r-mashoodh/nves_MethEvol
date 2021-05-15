## Functional annotation based on Blast2GO output

library(topGO)
library(tidyverse)

# ## use dm blast annotation
# 
# dm <- read.table("beetle_GO/GO_fly_clean.tab", sep="\t", quote=F)
# 
# dm.go <- dm %>% 
#   select(V1, V4) %>% 
#   rename(vesp.prot = V1, 
#          GO = V4)
# 
# ## need to get gene to prot mappings (might be in ortholog folder)
# 
# # contains protein ids (XP_)  
# proteins <- read.table("proteins_new.txt", stringsAsFactors = F) %>% 
#   distinct() %>% 
#   select(V1, V3)
# 
# # contains gene names attached to protein
# genes <- read.table("genes_new.txt", stringsAsFactors = F) %>% 
#   select(V1, V2)
# 
# genes.prot <- left_join(genes, proteins)
# colnames(genes.prot) <- c("LOC", "vesp.gene", "vesp.prot")
# # colnames(genes.prot) <- c("region_ID", "vesp.gene", "vesp.prot")
# # write.table(genes.prot, "../../Gen30_WGBS/Evol_Methylation/genes_to_protein_map.txt",
# #             quote=F, row.names=F, sep="\t")
# 
# # merge together to find a complete file
# dm.gene.prot <- left_join(dm.go, genes.prot) %>% 
#   distinct(vesp.gene, LOC, .keep_all = T)
# 
# # Create a topGO formatted annotation file, write it to csv for future
# dm.blast.GO <- dm.gene.prot %>% 
#   select(vesp.gene, GO)
# 
# write.table(dm.blast.GO, "orthologs/dm_blast2go.tab", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

### Gene lists ###
### These contain genes that are p<0.05 & abs(Log2FoldChange >= 1)
Evol.2fold <- read.csv("Evol_2fold_deseq2.csv")
FCenv.2fold <- read.csv("FCenv_2fold_deseq.csv")


### separating lists into lost vs gained
#### look at lost vs gained 
overlap.2fold <- intersect(Evol.2fold$gene, FCenv.2fold$gene)
gained.2fold <- setdiff(Evol.2fold$gene, FCenv.2fold$gene)
lost.2fold <- setdiff(FCenv.2fold$gene, Evol.2fold$gene)

##### TOPGO #####

## use topGO
geneID2GOdm <- readMappings(file="dm_blast2go.tab")
geneUniversedm <- names(geneID2GOdm)
EvolgenesOfInterest <- as.character(Evol.2fold$gene)
FCgenesOfInterest <- as.character(FCenv.2fold$gene)

### test of where difference are coming

evol.go <- left_join(Evol.2fold, dm.blast.GO, by = c("gene" = "vesp.gene"))
FCenv.go <- left_join(FCenv.2fold, dm.blast.GO, by = c("gene" = "vesp.gene"))


# identify where your genes of interest are within the gene Universe
EvolgeneList <- factor(as.integer(geneUniversedm %in% EvolgenesOfInterest))
names(EvolgeneList) <- geneUniversedm

FCgeneList <- factor(as.integer(geneUniversedm %in% FCgenesOfInterest))
names(FCgeneList) <- geneUniversedm

### Evol genes ###
Evol.myGOdata.dm <- new("topGOdata", 
                        description="Evol2fold", 
                        ontology="BP", 
                        allGenes=EvolgeneList,
                        annot = annFUN.gene2GO,
                        gene2GO = geneID2GOdm,
                        nodeSize=10)

Evol.resultFisher.dm <- runTest(Evol.myGOdata.dm, algorithm="weight01", statistic="fisher")


Evol.allRes.dm <- GenTable(Evol.myGOdata.dm, 
                           classicFisher = Evol.resultFisher.dm,
                           ranksOf = "classicFisher", 
                           topNodes = 200)

write.csv(Evol.allRes.dm, "Evol_GO_flyannot.csv")

FC.myGOdata <- new("topGOdata", 
                   description="FCenv2fold", 
                   ontology="BP", 
                   allGenes=FCgeneList,
                   annot = annFUN.gene2GO,
                   gene2GO = geneID2GOdm,
                   nodeSize=10)

FC.resultFisher <- runTest(FC.myGOdata, algorithm="classic", statistic="fisher")
FC.resultKS <- runTest(FC.myGOdata, algorithm="classic", statistic="ks")


FC.allRes <- GenTable(FC.myGOdata, 
                      classicFisher = FC.resultFisher,
                      ranksOf = "classicFisher", 
                      topNodes = 200)

write.csv(FC.allRes, "FCenv_GO_flyannot.csv")


#### look at lost vs gained 

overlap.2fold <- intersect(Evol.2fold$gene, FCenv.2fold$gene)
gained.2fold <- setdiff(Evol.2fold$gene, FCenv.2fold$gene)
lost.2fold <- setdiff(FCenv.2fold$gene, Evol.2fold$gene)


### LOST GENES ###
# identify where your genes of interest are within the gene Universe
LostList <- factor(as.integer(geneUniversedm %in% lost.2fold))
names(LostList) <- geneUniversedm

Lost.myGOdata <- new("topGOdata", 
                        description="Evol2fold", 
                        ontology="BP", 
                        allGenes=LostList,
                        annot = annFUN.gene2GO,
                        gene2GO = geneID2GOdm,
                     nodeSize=10)

Lost.resultFisher <- runTest(Lost.myGOdata, algorithm="weight01", statistic="fisher")

Lost.allRes <- GenTable(Lost.myGOdata, 
                           classicFisher = Lost.resultFisher,
                           ranksOf = "classicFisher", 
                           topNodes = 100)
show(Lost.allRes)
write.csv(Lost.allRes, "Lost_GO_flyannot.csv")


### OVERLAP ###

OverlapList <- factor(as.integer(geneUniversedm %in% overlap.2fold))
names(OverlapList) <- geneUniversedm

Overlap.myGOdata.dm <- new("topGOdata", 
                        description="Evol2fold",
                        ontology="BP",
                        allGenes=OverlapList,
                        annot = annFUN.gene2GO,
                        gene2GO = geneID2GOdm,
                        nodeSize=10)

Overlap.resultFisher.dm <- runTest(Overlap.myGOdata.dm, algorithm="weight01", statistic="fisher")

Overlap.allRes.dm <- GenTable(Overlap.myGOdata.dm, 
                           classicFisher = Overlap.resultFisher.dm,
                           ranksOf = "classicFisher", 
                           topNodes = 100)
show(Overlap.allRes.dm)
write.csv(Overlap.allRes.dm, "Overlap_GO_flyannot.csv")


### gained ####

GainedList <- factor(as.integer(geneUniversedm %in% gained.2fold))
names(GainedList) <- geneUniversedm

Gained.myGOdata.dm <- new("topGOdata", 
                           description="Evol2fold",
                           ontology="BP",
                           allGenes=GainedList,
                           annot = annFUN.gene2GO,
                           gene2GO = geneID2GOdm,
                          nodeSize=10)

Gained.resultFisher.dm <- runTest(Gained.myGOdata.dm, algorithm="weight01", statistic="fisher")

Gained.allRes.dm <- GenTable(Gained.myGOdata.dm, 
                              classicFisher = Gained.resultFisher.dm,
                              ranksOf = "classicFisher", 
                              topNodes = 100)
show(Gained.allRes.dm)
write.csv(Gained.allRes.dm, "Gained_GO_flyannot.csv")



## use topGO with "all_go" (combined annotation)
# geneID2GOdm <- readMappings("all_go_combined.tab")
# geneUniversedm <- names(geneID2GOdm)

