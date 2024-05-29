### Sarkies et al., (2024) Gene body methylation evolves during the sustained loss of parental care in the burying beetle.

#### 1. BS_processing - contains scripts for trimming, mapping and quantification of methylation data.
   - 1_trimming.sh - * *quality and adaptor trimming* *
   - 2_genome_setup.sh - * *genome prep* *
   - 3_bismark_align.sh - * *read alignment* *
   - 4_quant_meth.sh - * *quantification of methylation* *
    
#### 2. RNA_processing - contains scripts for trimpping, mapping and quantification of RNA sequencing data.
   - 1_setup.sh - * *software for hisat2 pipeline* *
   - 2_index.sh - * *genome indexing step* *
   - 3_processing.sh - * *read alignment* *
   - 4_fastqc.sh - * *fastqc* *
   
#### 3. R_scripts - contains scripts for analysing methylation and transcriptomic data.
   - 1_RNA_xpsn_Deseq2.R - * *RNA-seq analysis* *
   - 2_Cpg_evol.R - * *Differential mCpG analysis* *
   - 3_Mixtools_classification_DMGs.R - * *Differential gene methylation analysis* *
   - 4_GeneVsMe_v2.R - * *Methylated gene classification and association with gene expression (levels and variability) analysis* *
   - 5_CpG_correlation_genes.R - * *Quantification of methylation correlation amongst genes* *


