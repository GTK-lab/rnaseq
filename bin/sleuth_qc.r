# =================================================================== #
# LOAD LIBRARIES
# =================================================================== #
library(optparse)
library(rhdf5)
library(sleuth)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)

# =================================================================== #
# PARSE COMMAND-LINE PARAMETERS
# =================================================================== #
option_list <- list(
    make_option(c("-i", "--h5_dir"        ), type="character", default='kallisto', metavar="string" , help="the directory storing abundance.h5  from Kallisto."                                     ),
    make_option(c("-t", "--tx2gene"       ), type="character", default=NULL      , metavar="path"   , help="the corresponding Ensembl transcript and gene IDs."                                     ),
    make_option(c("-r", "--sample_suffix" ), type="character", default=''        , metavar="string" , help="Suffix to remove after sample name in columns e.g. '.rmDup.bam' if 'DRUG_R1.rmDup.bam'."),
    make_option(c("-o", "--outdir"        ), type="character", default='./'      , metavar="path"   , help="Output directory."                                                                      ),
    make_option(c("-p", "--outprefix"     ), type="character", default='sleuth'  , metavar="string" , help="Output prefix."                                                                         ),
    make_option(c("-c", "--cores"         ), type="integer"  , default=1         , metavar="integer", help="Number of cores."                                                                       )
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$h5_dir)){
    print_help(opt_parser)
    stop("Please provide a directory storing abundance.h5.", call.=FALSE)
}

# =================================================================== #
# PREPARE PARAMETERS: metadata; formula; tx2gene map
# =================================================================== #
# Paths of bundance.h5 files
h5_dir = list.files(opt$h5_dir, pattern = "abundance.h5", recursive = T, full.names = T)

# Prepare metadata
sample.name     = basename(dirname(h5_dir))
name_components = strsplit(sample.name, "_")
n_components    = length(name_components[[1]])
decompose       = n_components!=1 && all(sapply(name_components, length)==n_components)
metadata        = data.frame(sample.name, sample = sample.name, row.names=1)
if (decompose) {
    groupings        = as.data.frame(lapply(1:n_components, function(i) sapply(name_components, "[[", i)))
    names(groupings) = paste0("Group", 1:n_components)
    n_distinct       = sapply(groupings, function(grp) length(unique(grp)))
    groupings        = groupings[n_distinct!=1 & n_distinct!=length(sample.name)]
    if (ncol(groupings)!=0) {
        metadata <- cbind(metadata, groupings)
    } else {
        decompose <- FALSE
    }
}
metadata$path = h5_dir

# Read in tx2gene
tx2gene = read.delim(file = opt$tx2gene, header = FALSE, row.names = NULL, sep = "\t")
colnames(tx2gene) = c("target_id", "ens_gene", "ext_gene")

# Design formula
design = ~ 1

# =================================================================== #
# FIT THE SLEUTH MODEL
# =================================================================== #
so = sleuth_prep(metadata,
    full_model = design,
    target_mapping = tx2gene,
    read_bootstrap_tpm = TRUE,
    extra_bootstrap_summary = TRUE,
    transformation_function = function(x) log2(x + 0.5),
    num_cores = opt$cores
)
so = sleuth_fit(so)
oe = sleuth_wt(so, which_beta = '(Intercept)')
oe_results = sleuth_results(oe, test ='(Intercept)' )

results.name = paste0(opt$outprefix, ".tsv")
write.table(oe_results, file = results.name, sep = "\t", row.names = FALSE)

# =================================================================== #
# PLOT QC
# =================================================================== #
norm_counts = sleuth_to_matrix(oe, which_df = "obs_norm", which_units = "est_counts")
allgenes = nrow(norm_counts)

PlotFile <- paste(opt$outprefix,".plots.pdf",sep="")

pdf(file=PlotFile, onefile=TRUE, width=7, height=7)
# PCA
ntop = c(500, allgenes)
for (top_genes in ntop) {
    pca.data = prcomp(t(norm_counts[1:top_genes, ]))
  
    pca.plot = data.frame(pca.data$x, metadata)
    pca.var  = round(pca.data$sdev^2/sum(pca.data$sdev^2)*100)
  
    plot.subtitle = ifelse(top_genes == "500", paste("Top", top_genes, "Genes"), "All Genes")
  
     pl = ggplot(pca.plot, aes(PC1, PC2, label = paste(" ", sample, " "))) +
        geom_point() +
        geom_text(check_overlap = TRUE, vjust = 0.5, hjust = "inward") +
        xlab(paste0("PC1: ", round(pca.var[1]), "% variance")) +
        ylab(paste0("PC2: ", round(pca.var[2]), "% variance")) +
        labs(title = paste0("First PCs on data"), subtitle = plot.subtitle) +
        theme(legend.position = "top",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(color = "black", fill = NA, size = 1))
    print(pl)
  
    if (decompose) {
        pc_names = colnames(pca.data$x)
        long_pc <- reshape(pca.plot, varying=pc_names, direction="long", sep="", timevar="component", idvar="pcrow")
        long_pc <- subset(long_pc, component<=5)
        long_pc_grp <- reshape(long_pc, varying=names(groupings), direction="long", sep="", timevar="grouper")
        long_pc_grp <- subset(long_pc_grp, grouper<=5)
        long_pc_grp$component <- paste("PC", long_pc_grp$component)
        long_pc_grp$grouper <- paste0(long_pc_grp$grouper, c("st","nd","rd","th","th")[long_pc_grp$grouper], " prefix")
        pl <- ggplot(long_pc_grp, aes(x=Group, y=PC)) +
            geom_point() +
            stat_summary(fun=mean, geom="line", aes(group = 1)) +
            labs(x=NULL, y=NULL, subtitle = plot.subtitle, title="PCs split by sample-name prefixes") +
            facet_grid(component~grouper, scales="free_x") +
            scale_x_discrete(guide = guide_axis(n.dodge = 3))
        print(pl)
    }
}

# WRITE PC1 & PC2 VALUES TO FILE
pca.vals           = pca.data$x[, c("PC1", "PC2")]
colnames(pca.vals) = paste0(colnames(pca.vals), ":", pca.var[1:2], "% variance")
pca.vals           = cbind(sample = rownames(pca.vals), pca.vals)
write.table(pca.vals, file = paste0(opt$outprefix, ".pca.vals.txt"),
            row.names = FALSE, col.names = TRUE, sep = "\t", quote = TRUE)


# SAMPLE CORRELATION HEATMAP
sampleDists      = dist(t(norm_counts))
sampleDistMatrix = as.matrix(sampleDists)
colors           = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pheatmap(
  sampleDistMatrix,
  clustering_distance_rows=sampleDists,
  clustering_distance_cols=sampleDists,
  col=colors,
  main=paste("Euclidean distance samples")
)

write.table(cbind(sample = rownames(sampleDistMatrix), sampleDistMatrix),file=paste(opt$outprefix, ".sample.dists.txt", sep=""),
            row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

dev.off()

# =================================================================== #
# SAVE NORMALIZATION FACTORS
# =================================================================== #
norm.factors.dir = "norm_factors/"
if(file.exists(norm.factors.dir) == FALSE) {
    dir.create(norm.factors.dir, recursive = TRUE)
}
norm.factors.file = paste(norm.factors.dir, opt$outprefix, ".size_factor.RData", sep = "")

norm.factors = norm_factors(norm_counts)
save(norm.factors, file = norm.factors.file)

for(name in names(norm_factors(norm_counts))) {
    sizeFactorFile = paste0(norm.factors.dir, name, ".txt", sep = "")
    write(as.numeric(norm_factors(norm_counts)[name]), file = sizeFactorFile)
}


# =================================================================== #
# R SESSION INFO
# =================================================================== #
RLogFile = "R_sessionInfo.log"

sink(RLogFile)
a <- sessionInfo()
print(a)
sink()