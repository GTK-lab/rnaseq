# ===================================================================== #
#  Transcript-level estimates improve gene-level inferences
# ===================================================================== #
library(tximport)
library(SummarizedExperiment)

agrs = commandArgs(trailingOnly = TRUE)
# args: NULL, kallisto, kallisto_merge
if (length(args) < 2) {
    stop("Usage: kallisto_tximport.r <coldata> <salmon_out>", call. = FALSE)
}

coldata     = agrs[1] # NULL
path        = agrs[2] # kallisto
sample_name = args[3] # kallisto_merge

prefix  = sample_name
tx2gene = "kallisto_tx2gene.csv"
rowdata = read.csv(tx2gene) # colnames(rowdata): TXNAME, GENEID
tx2gene = rowdata

fns   = list.files(path, pattern = "abundance.tsv", recursive = T, full.names = T)
names = basename(dirname(fns)) # meta.id or prefix setted by customers
names(fns) = names

if (file.exists(coldata)) {
    coldata = read.csv(coldata, sep="\t")
    coldata = coldata[match(names, coldata[,1]),]
    coldata = cbind(files = fns, coldata)
} else {
    message("ColData not avaliable ", coldata)
    coldata = data.frame(files = fns, names = names)
}

txi = tximport (fns, type = "kallisto", txOut = TRUE)
rownames(coldata) = coldata[["names"]]
rownames(rowdata) = rowdata[["TXNAME"]]

# counts of transcripts
se = SummarizedExperiment(assays = list(
        counts = txi[["counts"]],
        abundance = txi[["abundance"]],
        length = txi[["length"]]),
    colData = DataFrame(coldata)
    rowData = data.frame(TX = rownames(txi[["counts"]]))
)

# counts of genes
if (!is.null(tx2gene)) {
    gi    = summarizeToGene(txi, tx2gene = tx2gene, ignoreTxVersion = TRUE)
    gi.ls = summarizeToGene(txi, tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)
    gi.s  = summarizeToGene(txi, tx2gene = tx2gene, countsFromAbundance = "scaledTPM", ignoreTxVersion = TRUE)
    growdata = data.frame(GENEID = unique(tx2gene[,2]))
    growdata = data.frame(GENEID = growdata[match(rownames(gi[[1]]), growdata[["GENEID"]]), ])

    gse   = SummarizedExperiment(assays = list(
            counts = gi[["counts"]],
            abundance = gi[["abundance"]],
            length = gi[["length"]]),
        colData = DataFrame(coldata),
        rowData = growdata)
    gse.ls = SummarizedExperiment(assays = list(
            counts = gi.ls[["counts"]],
            abundance = gi.ls[["abundance"]],
            length = gi.ls[["length"]]),
        colData = DataFrame(coldata),
        rowData = growdata)
    gse   = SummarizedExperiment(assays = list(
            counts = gi.s[["counts"]],
            abundance = gi.s[["abundance"]],
            length = gi.s[["length"]]),
        colData = DataFrame(coldata),
        rowData = growdata)
}
build_table = function(se.obj, slot) {
    cbind(rowData(se.obj)[, 1], assays(se.obj)[[slot]])
}

if(exists("gse")){
    write.table(build_table(gse, "abundance"), paste(c(prefix, "gene_tpm.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse, "counts"), paste(c(prefix, "gene_counts.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.ls, "abundance"), paste(c(prefix, "gene_tpm_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.ls, "counts"), paste(c(prefix, "gene_counts_length_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.s, "abundance"), paste(c(prefix, "gene_tpm_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
    write.table(build_table(gse.s, "counts"), paste(c(prefix, "gene_counts_scaled.tsv"), collapse="."), sep="\t", quote=FALSE, row.names = FALSE)
}

write.table(build_table(se,"abundance"), paste(c(prefix, "transcript_tpm.tsv"), collapse="."), sep="\t", quote=FALSE)
write.table(build_table(se, "counts"), paste(c(prefix, "transcript_counts.tsv"), collapse="."), sep="\t", quote=FALSE)