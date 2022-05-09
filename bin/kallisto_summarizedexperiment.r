library(SummarizedExperiment)

args = commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
    stop("Usage: kallisto_se.r <coldata> <counts> <tpm>", call.=FALSE)
}

coldata = args[1]
counts_fn = args[2]
tpm_fn = args[3]

tx2gene = "salmon_tx2gene.tsv"
info = file.info(tx2gene)
if (info$size == 0) {
    tx2gene = NULL
} else {
    rowdata = read.csv(tx2gene)
    tx2gene = rowdata
}

counts = read.csv(counts_fn, row.names=1, sep="\t")
tpm = read.csv(tpm_fn, row.names=1, sep="\t")

se = SummarizedExperiment(
    assays = list(
        counts = counts,
        abundance = tpm
    ),
    colData = DataFrame(coldata)
)
saveRDS(se, file = paste0(tools::file_path_sans_ext(counts_fn), "rds"))