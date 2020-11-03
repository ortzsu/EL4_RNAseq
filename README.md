# EL4_RNAseq
R script for RNAseq analysis using DESeq2

Three independent wild type “empty” and Furin super-enhancer deleted SE EL-4 cell clones were left unstimulated or stimulated with PMA/ionomycin for 6 hours. Total RNA was isolated and sent to NovoGene for RNA sequencing using the Illuminia PE150 platform.
Sequenced reads were mapped to the reference genome (GRCm38.p6.genome.fa) using STAR software (version 2.7.). Annotation file used: gencode.vM23.annotation.gtf
Raw counts were obtained using –quantMode GeneCounts option in STAR. The raw count matrix was imported into R and used as input for DESeq2 for normalization and differential gene expression analysis.
Genome_build: mm10
Data are deposited in GEO (accession number: GSE158456), .bed files are accessible at UCSC: https://genome.ucsc.edu/s/ortzsu/EL%2D4%20%22empty%22%20and%20SE%20clones
