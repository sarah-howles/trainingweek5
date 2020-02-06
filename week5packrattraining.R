getwd

#install tximport and DEseq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("rhdf5")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Mmusculus.v79")

library("tximport")
library("DESeq2")
library("rhdf5")
library("tidyverse")
library("biomaRt")

samples <- read_tsv("data/obds_sampletable.csv")
samples
transcripts_of_interest <- read_tsv("data/pseudoaligned/ERR1755082/abundance.tsv") 
transcripts_of_interest <- transcripts_of_interest$target_id

mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "useast.ensembl.org")
annot_genes <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id", "mgi_symbol"), 
                     values = transcripts_of_interest, mart = mouse)
 
gx <- EnsDb.Mmusculus.v79
gnm_tx <- transcripts(gx, return.type = "data.frame")
gnm_tx_filtered <- gnm_tx[c("tx_name", "gene_id")]


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

head(tx2gene)

files <- file.path(dir, "kallisto_boot", samples$run, "abundance.h5")
names(files) <- paste0("sample", 1:6)
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
