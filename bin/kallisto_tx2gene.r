# ===================================================================== #
#  Make Transcript to Gene Map
# ===================================================================== #

library(GenomicFeatures)

#---------------------------------------------------------------------- #
#  Create tx2gene object with GTF files
#---------------------------------------------------------------------- #
args = commandArgs(trailingOnly = TRUE)
gtf  = args
txdb = makeTxDbFromGFF(file = gtf, format = "gtf")
k    = keys(txdb, keytype = "TXNAME")
tx2gene = select(txdb, k, "GENEID", "TXNAME")
write.csv(tx2gene, "kallisto_tx2gene.csv", row.names = FALSE)

