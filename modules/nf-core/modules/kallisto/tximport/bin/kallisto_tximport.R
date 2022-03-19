#!/usr/bin/env Rscript

library(tximport)

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop("Usage: kallisto_tximport.R <tx2gene_file> <kallisto_dir>", call.=FALSE)
}

tx2geneFile = args[1]
file_path = args[2]

kallisto_filepaths=file.path(file_path,dir(file_path,recursive=TRUE, pattern="*_abundance.txt"))

# This hack below is needed to ensure the columns get ordered as
# Sample_1, Sample_2, Sample_3 ..... Sample_10, Sample_11, Sample_12
# Without this hack, the columns will be Sample_1, Sample_10, Sample_11, Sample_12, Sample_2, Sample_3

# For the above reasons, I have column 1 in my original samplesheet to always have a running number (_\d+) after sample_id
# sample_id,fastq_1,fastq_2
# Sample_1,<path>,<path>
# Sample_2,<path>,<path>
# Sample_3,<path>,<path>
# ......
# Sample_10,<path>,<path>
# Sample_11,<path>,<path>

kallisto_filepaths=kallisto_filepaths[order(as.numeric(gsub(".*?kallisto/.*?_(\\d+).kallisto_abundance.txt","\\1", kallisto_filepaths)))]

samples = data.frame(samples = gsub(".*?kallisto/(.*?).kallisto_abundance.txt", "\\1", kallisto_filepaths) )

row.names(samples)=samples[,1]

names(kallisto_filepaths) = samples$samples

my_tx2gene=read.csv(tx2geneFile,sep = "\t",stringsAsFactors = F, header=F)

kallisto_txi_default <- tximport(kallisto_filepaths, type = "kallisto", tx2gene = my_tx2gene, ignoreAfterBar = TRUE, countsFromAbundance="no")
kallisto_txi_lengthScaledTPM <- tximport(kallisto_filepaths, type = "kallisto", tx2gene = my_tx2gene, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")
kallisto_txi_scaledTPM <- tximport(kallisto_filepaths, type = "kallisto", tx2gene = my_tx2gene, ignoreAfterBar = TRUE, countsFromAbundance="scaledTPM")
kallisto_txi_dtuScaledTPM <- tximport(kallisto_filepaths, type = "kallisto", tx2gene = my_tx2gene, ignoreAfterBar = TRUE, countsFromAbundance="dtuScaledTPM", txOut = TRUE)

#Gene Count Summarization
write.csv(as.data.frame(kallisto_txi_default$counts),
                  file = "Kallisto_GeneCount_default.csv",
                  quote = FALSE)

write.csv(as.data.frame(kallisto_txi_lengthScaledTPM$counts),
                  file = "Kallisto_GeneCount_lengthScaledTPM.csv",
                  quote = FALSE)

write.csv(as.data.frame(kallisto_txi_scaledTPM$counts),
                  file = "Kallisto_GeneCount_scaledTPM.csv",
                  quote = FALSE)

write.csv(as.data.frame(kallisto_txi_scaledTPM$counts),
                  file = "Kallisto_GeneCount_scaledTPM.csv",
                  quote = FALSE)

write.csv(as.data.frame(kallisto_txi_dtuScaledTPM$counts),
                  file = "Kallisto_GeneCount_dtuScaledTPM.csv",
                  quote = FALSE)

# TPM

write.csv(as.data.frame(kallisto_txi_default$abundance),
                  file = "Kallisto_Abundance_default.csv",
                  quote = FALSE)

write.csv(as.data.frame(kallisto_txi_lengthScaledTPM$abundance),
                  file = "Kallisto_Abundance_lengthScaledTPM.csv",
                  quote = FALSE)

write.csv(as.data.frame(kallisto_txi_scaledTPM$abundance),
                  file = "Kallisto_Abundance_scaledTPM.csv",
                  quote = FALSE)

write.csv(as.data.frame(kallisto_txi_dtuScaledTPM$abundance),
                  file = "Kallisto_Abundance_dtuScaledTPM.csv",
                  quote = FALSE)
