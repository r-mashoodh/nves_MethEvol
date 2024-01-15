# Using mixture models to classify methylated vs unmethylated vespilloides genes
# https://tinyheero.github.io/2015/10/13/mixture-model.html
# https://jbhender.github.io/Stats506/F18/GP/Group7.html

library(tidyverse)
library(mixtools)
library(broom)
library(gridExtra)

## Read in raw data
weighted <- read.table("Nves.weighted.bedGraph", sep="\t", na.strings = c("."))

## Clean up data table and filter down to **genes only**
weighted.genes <- weighted %>%
  select(V1, V3, V4, V5, V9, V10, V11, V12) %>%
  rename(contig = V1,
         region_type = V3,
         start = V4,
         stop = V5,
         region_ID = V9,
         meth_avg = V10, #avg across the CpGs present in the region
         CpG_count = V11, #this is the # of CpGs
         sample_ID = V12) %>%
  mutate(group = substr(sample_ID, 1, 4),
         length = stop - start,
         evol = substr(sample_ID, 1, 2)) %>%
  filter(region_type == "gene")

## move to a wider format
weighted.wide <- weighted.genes %>%
  select(region_ID, sample_ID, meth_avg) %>%
  pivot_wider(names_from = sample_ID, values_from = meth_avg)

## look at NAs
weighted.wide$na_count <- apply(weighted.wide[,2:8], 1, function(x) sum(is.na(x)))
weighted.wide$na_count2 <- apply(weighted.wide[,9:15], 1, function(x) sum(is.na(x)))

## filter by number of NAs
wide.filtered <- weighted.wide %>%
  filter(na_count <= 1 & na_count2 <= 1) %>% #allows for no missing data across groups
  mutate(mean_meth = rowMeans(select(., 2:15), na.rm=T)) %>%
  select(region_ID, mean_meth) %>%
  column_to_rownames(var = "region_ID")

mat.wide <- data.matrix(wide.filtered)

weighted.wide %>%
  filter(na_count <= 1 & na_count2 <= 1) %>% #allows for no missing data across groups
  mutate(mean_meth = rowMeans(select(., 2:15), na.rm=T)) %>%
  select(region_ID, mean_meth) %>% 
  write_csv(., file="nves_all_mCpG.csv")

## use mixtools 
mixmdl <- normalmixEM(mat.wide)

## plot model
plot(mixmdl, which=2)
lines(density(mixmdl$x), lty=2, lwd=0.8)

mixmdl$mu
# [1] 0.4308975 9.6761851

mixmdl$sigma
# [1] 0.09499342 11.67491285

mixmdl$lambda #height
# [1] 0.6095893 0.3904107
# can also interpret this as mixing weights, 60% of genes are not methylated

# Create a df of probabilities
post.df <- as.data.frame(cbind(x = mixmdl$x, mixmdl$posterior))
post.df <- cbind(mat.wide, post.df) %>% 
  rownames_to_column(var="region_ID") %>% 
  select(-x)

head(post.df)
# region_ID mean_meth       comp.1      comp.2
# 1 LOC108556273 0.5325093 9.934748e-01 0.006525172
# 2 LOC108557084 0.4329244 9.961859e-01 0.003814135
# 3 LOC108564750 0.3940229 9.959139e-01 0.004086050
# 4 LOC108560513 0.9968625 1.512748e-05 0.999984873
# 5 LOC108562766 0.3257794 9.932727e-01 0.006727338
# 6 LOC108561911 1.7680168 1.116265e-38 1.000000000

# Use ther posterior probs to check # of genes in each category
length(which(post.df$comp.1 < 0.01))
#[1] 4431
# this is the number of genes *not* in comp1 (i.e. methylated)

length(which(post.df$comp.2 < 0.01))
#[1] 5721
# this is the num of genes *not* in comp2 (i.e. not methylated)

nrow(mat.wide)
#[1] 11048

### Plot genes at different cutoffs
## need to go a lot lower that 0.05 to start to see differences in the number of methylated gens
post.df %>%
  mutate(label = ifelse(comp.1 < 0.01, 1, 0)) %>% 
  ggplot(aes(x = factor(label))) +
  geom_bar() +
  xlab("Component") +
  ylab("Number of Data Points")

## 
summary(post.df$comp.1)
summary(post.df$comp.2)

## plot density
post.df %>%
  mutate(label = ifelse(comp.1 < 0.01, 1, 0)) %>% 
  ggplot(aes(x = mean_meth, group=as.factor(label), fill=as.factor(label))) +
  geom_density(alpha=0.5, adjust=2)

#### plot mixtools

#### now run a test to see which are differentially methylated
### using binomial

# get data that contains meth and unmeth counts
nves.meth.all <- read.table("Nves.allsamples.bedGraph", sep="\t", na.strings = c(".")) 

## filter methylated genes
methylated.mixmdl <- post.df %>% 
  filter(comp.1 < 0.01)

qplot(methylated.mixmdl$mean_meth, geom="histogram")

# cleanup dataframe and filter by 'methylated' genes
counts.for.mixmdl <- nves.meth.all %>% 
  select(V1, V3, V4, V5, V9, V11, V12, V13) %>% 
  dplyr::rename(contig = V1,
         region_type = V3,
         start = V4, 
         stop = V5,
         region_ID = V9,
         meth_count = V11,
         unmeth_count = V12,
         sample_ID = V13) %>% 
  filter(region_type=="gene") %>%
  mutate(group = substr(sample_ID, 1, 4),
         length = stop - start,
         evol = substr(sample_ID, 1, 2),
         env = substr(sample_ID, 3, 4),
         evol_num = ifelse(evol=="FC", 0, 1),
         env_num = ifelse(env == "FC", 0, 1),
         meth_mean = meth_count / (meth_count + unmeth_count)) %>%
  filter(region_ID %in% methylated.mixmdl$region_ID)

# # including all samples plus an interaction (factor coded)
# binom.mixmdl <- counts.for.mixmdl %>% 
#   group_by(region_ID) %>% 
#   do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$evol*.$env, family="binomial"))) %>%
#   filter(term == ".$evolNC")
# 
# binom.mixmdl$padj <- p.adjust(binom.mixmdl$p.value, n=nrow(binom.mixmdl), method="BH")
# 
# hist(binom.mixmdl$p.value)
# 
# length(which(binom.mixmdl$padj < 0.01))
# # [1] 1429
# 
# length(which(binom.mixmdl$padj < 0.05))
# # [1] 1891


# ### add env as a covariate?
# 
# cov.test <- counts.for.mixmdl %>% 
#   group_by(region_ID) %>% 
#   do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$evol_num + .$env_num, family="binomial"))) %>% 
#   filter(term == ".$evol_num")
# 
# cov.test$padj <- p.adjust(cov.test$p.value, method="BH", n=nrow(cov.test))
# 
# cov.filtered <- cov.test %>% 
#   filter(padj < 0.01) %>% 
#   arrange(-abs(estimate))
# 
# # this is actually an interaction
# counts.for.mixmdl %>% 
#   filter(region_ID == "LOC108563366") %>% 
#   ggplot(., aes(x=group, y=meth_mean)) +
#   geom_boxplot()
# 
# # barely detectable effect?
# counts.for.mixmdl %>% 
#   filter(region_ID == "LOC108558176") %>% 
#   ggplot(., aes(x=group, y=meth_mean)) +
#   geom_boxplot()
# 
# ## so really adding env as a covariate doesn't really work 
# ## need to just do it separately, so why not filter by group

binom.mixmdl.evol <- counts.for.mixmdl %>% 
  filter(group == "FCFC" | group == "NCNC") %>% 
  group_by(region_ID) %>% 
  do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$evol, family="binomial"))) %>%
  filter(term == ".$evolNC")

binom.mixmdl.evol$padj <- p.adjust(binom.mixmdl.evol$p.value, 
                                   n=nrow(binom.mixmdl.evol), method="BH")

hist(binom.mixmdl.evol$p.value)

length(which(binom.mixmdl.evol$padj < 0.01))
# [1] 1530
length(which(binom.mixmdl.evol$padj < 0.05))

# now try to incorporate mean diff
group_means <- counts.for.mixmdl %>% 
  group_by(region_ID, group) %>% 
  summarise(gene_mean = mean(meth_mean, na.rm=T)) %>% 
  pivot_wider(names_from = group, values_from = gene_mean)

## make a plot of % diff and p-value (-log10(p))
## volcano
binom.mixmdl.evol %>% 
  left_join(group_means) %>% 
  mutate(FCFC_NCNC = (NCNC - FCFC)*100) %>% 
  ggplot(aes(x=FCFC_NCNC, y=-log10(padj), color=(abs(FCFC_NCNC) > 15 & -log10(padj) > 2))) +
  geom_point() +
  geom_vline(xintercept=c(-15, +15), col="blue", linetype="dotted") +
  geom_hline(yintercept=2, col="blue", linetype="dotted") +
  scale_color_manual(values=c("grey", "black")) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text = element_text(size=20))
  ggtitle("evol.filtered volcano")
  
# get a table of sig p-values
evol.filtered <- binom.mixmdl.evol %>% 
  filter(padj < 0.01) %>% 
  arrange(-abs(estimate)) %>% 
  left_join(group_means) %>% 
  mutate(FCFC_NCNC = (NCNC - FCFC)*100) %>% 
  filter(abs(FCFC_NCNC) >= 15) ## > 15% difference

evol.filtered.nothresh <- binom.mixmdl.evol %>% 
  filter(padj < 0.01)

nrow(evol.filtered)
# [1] 161

test.evol <- counts.for.mixmdl %>% 
  group_by(region_ID) %>% 
  do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$evol, family="binomial"))) %>%
  filter(term == ".$evolNC")

test.evol$padj <- p.adjust(test.evol$p.value, method="BH", n=nrow(test.evol))
length(which(test.evol$padj < 0.01))

###### permutation for pop evol comparison ####
# permute_evol <- map_dfr(1:1000, ~ counts.for.mixmdl %>%
#                           group_by(region_ID) %>%
#                           mutate(evol_num = evol_num[sample(row_number())]) %>% # permutation shuffle the 'group' column within each 'level'.
#                           do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$evol_num, family="binomial"))) %>%
#                           filter(term == ".$evol_num") %>%
#                           mutate(padj = p.adjust(p.value, method="BH", n=nrow(binom.mixmdl.evol))) %>%
#                           ungroup() %>%
#                           summarise(p = sum(padj < 0.01)))
# 
# 
# hist(permute_evol$p)
# summary(permute_evol$p)
# 

#write.csv(permute_evol, "1000_permutations_mixmdl_binomial.csv")

## this might have been created with a subsetted df
permute1 <- read.csv("1000_permutations_mixmdl_binomial_v1.csv")

hist(permute1$p)
summary(permute1$p)

permute1 %>%
  ggplot(., aes(x=p)) +
  geom_histogram(bins=75) +
  geom_vline(xintercept = 1486, color="red") +
  geom_vline(xintercept = 1515, color="blue") +
  theme_bw()


# ## how many overlap with the plasticity contrasts
# length(intersect(evol.filtered$region_ID, env.filtered$region_ID))
# # 690 overlaps with FCplasticity
# 
# length(intersect(evol.filtered$region_ID, nc.env.filtered$region_ID))
# # 862 overlaps with NC plasticity
# 
# # these are the genes that are consistently diff methylated in a FC environment
# plasticity.meth <- intersect(nc.env.filtered$region_ID, env.filtered$region_ID)
# 
# # remove plasticity genes in either 
# evol.genes.no.env <- evol.filtered %>% 
#   filter(!region_ID %in% nc.env.filtered$region_ID) %>% 
#   filter(!region_ID %in% env.filtered$region_ID)
# # this leaves 322 lets plot and see what they look like
# 
# evol.only.genes <- evol.genes.no.env$region_ID[1:9]
# 
# plot_gene <- function(gene = NULL) {
#   counts.for.mixmdl %>% 
#     filter(group == "FCFC" | group == "NCNC") %>% 
#     filter(region_ID == gene) %>% 
#     ggplot(., aes(x=group, y=meth_mean, color=evol)) +
#     geom_boxplot() +
#     theme_bw() +
#     ggtitle(gene)
# }
# 
# grid.arrange(
#   plot_gene(evol.only.genes[1]),
#   plot_gene(evol.only.genes[2]),
#   plot_gene(evol.only.genes[3]),
#   plot_gene(evol.only.genes[4]),
#   plot_gene(evol.only.genes[5]),
#   plot_gene(evol.only.genes[6]),
#   plot_gene(evol.only.genes[7]),
#   plot_gene(evol.only.genes[8]),
#   ncol=3
# )


####### test for env (within FC-pop)

binom.mixmdl.env <- counts.for.mixmdl %>% 
  filter(group == "FCFC" | group == "FCNC") %>% 
  group_by(region_ID) %>% 
  do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$env_num, family="binomial"))) %>%
  filter(term == ".$env_num")

binom.mixmdl.env$padj <- p.adjust(binom.mixmdl.env$p.value, n=nrow(binom.mixmdl.env), method="BH")

hist(binom.mixmdl.env$p.value)

length(which(binom.mixmdl.env$padj < 0.01))
# [1] 1486

env.filtered <- binom.mixmdl.env %>% 
  filter(padj < 0.01) %>% 
  arrange(-abs(estimate)) %>% 
  left_join(group_means) %>% 
  mutate(FCFC_FCNC = (FCNC - FCFC)*100) %>% 
  filter(abs(FCFC_FCNC) > 15) ## > 15% difference

env.filtered.nothresh <- binom.mixmdl.env %>% 
  filter(padj < 0.01)

binom.mixmdl.env %>%
  left_join(group_means) %>% 
  mutate(FCFC_FCNC = (FCNC - FCFC)*100) %>% 
  ggplot(aes(x=FCFC_FCNC, y=-log10(padj), 
             color=(abs(FCFC_FCNC) > 15 & -log10(padj) > 2))) +
  geom_point() +
  geom_vline(xintercept=c(-15, 15), col="blue", linetype="dotted") +
  geom_hline(yintercept=2, col="blue", linetype="dotted") +
  scale_color_manual(values=c("grey", "black")) +
  theme(legend.position = "none") +
  ggtitle("env.filtered volcano")

nrow(env.filtered)
# 160

# whats overlapping
length(intersect(evol.filtered$region_ID, env.filtered$region_ID))
# [1] 50

# are they in the same direction?
# what's new?


####### test for env (within NC-pop)

binom.mixmdl.env.nc <- counts.for.mixmdl %>% 
  filter(group == "NCFC" | group == "NCNC") %>% 
  group_by(region_ID) %>% 
  do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$env_num, family="binomial"))) %>%
  filter(term == ".$env_num")

binom.mixmdl.env.nc$padj <- p.adjust(binom.mixmdl.env.nc$p.value, n=nrow(binom.mixmdl.env.nc), method="BH")

hist(binom.mixmdl.env.nc$p.value)

length(which(binom.mixmdl.env.nc$padj < 0.01))
# [1] 1788

nc.env.filtered <- binom.mixmdl.env.nc %>% 
  filter(padj < 0.01) %>%
  left_join(group_means) %>% 
  mutate(NCFC_NCNC = (NCNC - NCFC)*100) %>% 
  filter(abs(NCFC_NCNC) > 15) %>% ## > 15% difference
  arrange(-abs(NCFC_NCNC)) 

nc.env.filtered.nothresh <- binom.mixmdl.env.nc %>% 
  filter(padj < 0.01)

nrow(nc.env.filtered)
#207

binom.mixmdl.env.nc %>%
  left_join(group_means) %>% 
  mutate(NCFC_NCNC = (NCNC - NCFC)*100) %>% 
  ggplot(aes(x=NCFC_NCNC, y=-log10(padj), 
             color=(abs(NCFC_NCNC) > 15 & -log10(padj) > 2))) +
  geom_point() +
  geom_vline(xintercept=c(-15, 15), col="blue", linetype="dotted") +
  geom_hline(yintercept=2, col="blue", linetype="dotted") +
  scale_color_manual(values=c("grey", "black")) +
  theme(legend.position = "none") +
  ggtitle("NCenv.filtered volcano")


# whats overlapping with FC env in FC-pop
length(intersect(nc.env.filtered$region_ID, env.filtered$region_ID))
# [1] 40

# are they in the same direction?

# what's not overlapping
length(setdiff(evol.filtered$region_ID, env.filtered$region_ID))
# [1] 96


# #### FCFC vs NCFC 
# 
# binom.FCFC_NCFC <- counts.for.mixmdl %>% 
#   filter(group == "FCFC" | group == "NCFC") %>% 
#   group_by(region_ID) %>% 
#   do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$evol_num, family="binomial"))) %>%
#   filter(term == ".$evol_num")
# 
# binom.FCFC_NCFC$padj <- p.adjust(binom.FCFC_NCFC$p.value, 
#                                  n=nrow(binom.FCFC_NCFC), method="BH")
# 
# hist(binom.FCFC_NCFC$p.value)
# 
# FCFC_NCFC.filtered <- binom.FCFC_NCFC %>% 
#   filter(padj < 0.01) %>% 
#   arrange(-abs(estimate)) %>% 
#   left_join(group_means) %>% 
#   mutate(FCFC_NCFC = (NCFC - FCFC)*100) %>% 
#   filter(abs(FCFC_NCFC) > 15) ## > 15% difference
# 
# #####
# 
# binom.FCNC_NCNC <- counts.for.mixmdl %>% 
#   filter(group == "FCFC" | group == "NCNC") %>% 
#   group_by(region_ID) %>% 
#   do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$evol_num, family="binomial"))) %>%
#   filter(term == ".$evol_num")
# 
# binom.FCNC_NCNC$padj <- p.adjust(binom.FCNC_NCNC$p.value, 
#                                  n=nrow(binom.FCNC_NCNC), method="BH")
# 
# hist(binom.FCNC_NCNC$p.value)
# 
# FCNC_NCNC.filtered <- binom.FCNC_NCNC %>% 
#   filter(padj < 0.01) %>% 
#   arrange(-abs(estimate)) %>% 
#   left_join(group_means) %>% 
#   mutate(FCNC_NCNC = (NCNC - FCNC)*100) #%>% 
#   #filter(abs(FCNC_NCNC) > 15) ## > 15% difference
# 
# ## assuming main effects are the overlap 
# main.effects <- intersect(FCFC_NCFC.filtered$region_ID, FCNC_NCNC.filtered$region_ID)
# 
# plot_gene <- function(gene = NULL) {
#   counts.for.mixmdl %>%
#     #filter(group == "FCFC" | group == "NCNC") %>%
#     filter(region_ID == gene) %>%
#     ggplot(., aes(x=group, y=meth_mean, color=evol)) +
#     geom_boxplot() +
#     theme_bw() +
#     ggtitle(gene)
# }
# 
# grid.arrange(
#   plot_gene(main.effects[1]),
#   plot_gene(main.effects[2]),
#   plot_gene(main.effects[3]),
#   plot_gene(main.effects[4]),
#   plot_gene(main.effects[5]),
#   plot_gene(main.effects[6]),
#   ncol=3
# )
# 
# # What genes have no effect of env in either pop?
# # get genes that have no efffect padj > 0.10
# no.effect.FCenv <- binom.mixmdl.env %>% 
#   filter(padj > 0.1) %>% 
#   filter(region_ID %in% main.effects)
# 
# no.effect.NCenv <- binom.mixmdl.env.nc %>% 
#   filter(padj > 0.1) %>% 
#   filter(region_ID %in% main.effects) 
# 
# 
# grid.arrange(
#   plot_gene(all$region_ID[1]),
#   plot_gene(all$region_ID[2]),
#   plot_gene(all$region_ID[3]),
#   plot_gene(all$region_ID[4]),
#   ncol=2
# )
# 
# #############################################################################
# #test covariate or interaction
# test.full <- counts.for.mixmdl %>%
#   group_by(region_ID) %>%
#   do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$evol_num*.$env_num, family="binomial")))
# 
# n_genes <- length(unique(test.full$region_ID))
# 
# full.results <- test.full %>%
#   select(p.value, term, region_ID) %>%
#   pivot_wider(names_from = term, values_from=p.value) %>%
#   mutate(padj_Evol = p.adjust(`.$evol_num`, method="BH", n=n_genes),
#          padj_Env = p.adjust(`.$env_num`, method="BH", n=n_genes), 
#          padj_Int = p.adjust(`.$evol_num:.$env_num`, method="BH", n=n_genes)) %>% 
#   mutate(main_Evol = ifelse(padj_Evol < 0.01, 1, 0),
#          main_Env = ifelse(padj_Env < 0.01, 1, 0),
#          main_Int = ifelse(padj_Int < 0.01, 1, 0))
# 
# base_means <- wide.filtered %>% 
#   rownames_to_column(var = "region_ID")
# 
# evol_means <- counts.for.mixmdl %>% 
#   group_by(region_ID, evol) %>% 
#   summarise(gene_mean = mean(meth_mean)) %>% 
#   pivot_wider(names_from = evol, values_from = gene_mean) %>% 
#   mutate(diff = NC - FC)
# 
# Evolved_only <- full.results %>% 
#   filter(main_Int == 0 & main_Env == 0 & main_Evol == 1) %>%
#   left_join(evol_means) %>% 
#   arrange(-abs(diff))
# 
# plot_gene <- function(gene = NULL) {
#   counts.for.mixmdl %>%
#     #filter(group == "FCFC" | group == "NCNC") %>%
#     filter(region_ID == gene) %>%
#     ggplot(., aes(x=evol, y=meth_mean, color=evol)) +
#     geom_boxplot() +
#     theme_bw() +
#     ggtitle(gene)
# }
# 
# 
# evolved <- Evolved_only$region_ID[1:8]
# 
# grid.arrange(
#   plot_gene(evolved[1]),
#   plot_gene(evolved[2]),
#   plot_gene(evolved[3]),
#   plot_gene(evolved[4]),
#   plot_gene(evolved[5]),
#   plot_gene(evolved[6]),
#   plot_gene(evolved[7]),
#   plot_gene(evolved[8]),
#   ncol=4
# )
# 
# 
# 
# 
# ####################################################################################
# ### the key here is if FCFC is different from NCFC ... ie. do they become more similar once 
# ## care is restored?
# ## these are genes that behave similarly in either population
# plasticity <- intersect(env.filtered$region_ID, nc.env.filtered$region_ID)[1:6]
# 
# grid.arrange(
#   plot_gene(plasticity[1]),
#   plot_gene(plasticity[2]),
#   plot_gene(plasticity[3]),
#   plot_gene(plasticity[4]),
#   plot_gene(plasticity[5]),
#   plot_gene(plasticity[6]),
#   ncol=3
# )
# 
# 
# ## these are genes that behave similarly in either population
# nc.plasticity <- setdiff(nc.env.filtered$region_ID, env.filtered$region_ID)[1:6]
# 
# grid.arrange(
#   plot_gene(nc.plasticity[1]),
#   plot_gene(nc.plasticity[2]),
#   plot_gene(nc.plasticity[3]),
#   plot_gene(nc.plasticity[4]),
#   plot_gene(nc.plasticity[5]),
#   plot_gene(nc.plasticity[6]),
#   ncol=3
# )
# 
# 
# 
# binom.mixmdl.pop <- counts.for.mixmdl %>% 
#   filter(group == "FCFC" | group == "NCFC") %>% 
#   group_by(region_ID) %>% 
#   do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$evol_num, family="binomial"))) %>%
#   filter(term == ".$evol_num")
# 
# binom.mixmdl.pop$padj <- p.adjust(binom.mixmdl.pop$p.value, n=nrow(binom.mixmdl.pop), method="BH")
# 
# hist(binom.mixmdl.pop$p.value)
# 
# length(which(binom.mixmdl.pop$padj < 0.01))
# # [1] 1682
# 
# nc.env.pop <- binom.mixmdl.pop %>% 
#   filter(padj < 0.01) %>% 
#   arrange(-abs(estimate))
# 
# pop.in.FC <- nc.env.pop$region_ID[1:8]
# 
# grid.arrange(
#   plot_gene(pop.in.FC[1]),
#   plot_gene(pop.in.FC[2]),
#   plot_gene(pop.in.FC[3]),
#   plot_gene(pop.in.FC[4]),
#   plot_gene(pop.in.FC[5]),
#   plot_gene(pop.in.FC[6]),
#   plot_gene(pop.in.FC[7]),
#   plot_gene(pop.in.FC[8]),
#   ncol=4
# )
# 
# 
# ### FCNC vs NCNC
# binom.mixmdl.popNC <- counts.for.mixmdl %>% 
#   filter(group == "FCNC" | group == "NCNC") %>% 
#   group_by(region_ID) %>% 
#   do(tidy(glm(cbind(.$meth_count,.$unmeth_count) ~ .$evol_num, family="binomial"))) %>%
#   filter(term == ".$evol_num")
# 
# binom.mixmdl.popNC$padj <- p.adjust(binom.mixmdl.popNC$p.value, n=nrow(binom.mixmdl.popNC), method="BH")
# 
# hist(binom.mixmdl.popNC$p.value)
# 
# length(which(binom.mixmdl.popNC$padj < 0.01))
# # [1] 1682
# 
# nc.env.popNC <- binom.mixmdl.popNC %>% 
#   filter(padj < 0.01) %>% 
#   arrange(-abs(estimate))
# 
# pop.in.NC <- nc.env.popNC$region_ID[1:8]
# 
# grid.arrange(
#   plot_gene(pop.in.NC[1]),
#   plot_gene(pop.in.NC[2]),
#   plot_gene(pop.in.NC[3]),
#   plot_gene(pop.in.NC[4]),
#   plot_gene(pop.in.NC[5]),
#   plot_gene(pop.in.NC[6]),
#   plot_gene(pop.in.NC[7]),
#   plot_gene(pop.in.NC[8]),
#   ncol=4
# )
# 
# 
# ####### main effect genes
# 
# 
# 
# 
# #####################################################################################
# ################################### JUNK ############################################
# 
# 
# 
# #########
# a.gene <- counts.for.mixmdl %>% 
#   filter(region_ID == "LOC108562077")
# 
# mod <- glm(cbind(meth_count,unmeth_count) ~ evol_num*env, data=a.gene, family="binomial")
# summary(mod)
# 
# a.gene %>% 
#   ggplot(., aes(x=group, y=meth_mean)) +
#   geom_boxplot()
# 
# 
# counts.for.mixmdl %>% 
#   filter(region_ID == "LOC108564761") %>% 
#   ggplot(., aes(x=group, y=meth_mean)) +
#   geom_boxplot()
# 
# counts.for.mixmdl %>%
#   filter(region_ID == "LOC108569378") %>% 
#   ggplot(., aes(x=group, y=meth_mean)) +
#   geom_boxplot()
# 
# counts.for.mixmdl %>%
#   filter(region_ID == "LOC108556748") %>% 
#   ggplot(., aes(x=group, y=meth_mean)) +
#   geom_boxplot()
# 
# 
# env.term %>% 
#   filter(region_ID %notin% int.term$region_ID) %>% 
#   filter(abs(estimate) > 1)
# 
# 
# 
# 
# #permute
# 

# 
# ### try running ttest??
# 
# ttest.mixmdl <- counts.for.mixmdl %>% 
#   group_by(region_ID) %>% 
#   do(tidy(t.test(.$meth_mean ~ .$evol_num, conf.level=0.95)))
#      
# ttest.mixmdl$padj <- p.adjust(ttest.mixmdl$p.value, n=nrow(ttest.mixmdl), method="BH")
# 
# hist(ttest.mixmdl$p.value)
# 
# length(which(ttest.mixmdl$padj < 0.05))
# 




