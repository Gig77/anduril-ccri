library(optparse)

option_list <- list(
  make_option("--gsea-result-file", type="character", help="CSV file containing enrichment results produced by GSEA (e.g. 'gsea_report_for_na_pos_12345678.xls')."),
  make_option("--phenotype", type="character", help="either 'up' or 'down'"),
  make_option("--gene-set-dir", type="character", help="Directory containing individual gene set reports prodcued by GSEA."),
  make_option("--num-genes", type="integer", help="Number of genes in input .rnk file used for GEA analysis"),
  make_option("--categories", type="character", help="CSV file with category names for gene sets."),
  make_option("--output-file", type="character", help="Name of output file.")
)
opt <- parse_args(OptionParser(option_list=option_list))

# debug
# opt <- data.frame('gsea-result-file' = "/mnt/projects/iamp/results/anduril/execute/GSEA/case1/component/enrichedUp.csv", 'phenotype' = "up", 'gene-set-dir' = "/mnt/projects/iamp/results/anduril/execute/GSEA/case1/component/gsea_output", 'num-genes' = 87, 'categories' = '/mnt/projects/iamp/results/geneset_category.tsv', stringsAsFactors=F, check.names=F)   

stopifnot(opt$phenotype == 'up' || opt$phenotype == 'down')

print(paste("Adding core gene names with ranks to file", opt$'gsea-result-file', "..."))

# read GSEA results
gsea.result <- read.delim(opt$'gsea-result-file', stringsAsFactors = FALSE, check.names = FALSE)
out <- gsea.result[,!names(gsea.result) %in% c("")] # drop extra empty column
out$CORE_GENES <- NA

# read gene set to category assignments (if provided)
if (!is.null(opt$categories) && !is.na(opt$categories) && nchar(opt$categories) > 0) {
  categories <- read.delim(opt$categories, stringsAsFactors = FALSE)
  names(categories)[2] <- "CATEGORY"
  out <- merge(out, categories, by.x="NAME", by.y="geneset", all.x=T)
}


# for each significant gene set, load gene set result file and determine core enrichment gene names
processed <- 0
skipped <- 0
for (i in 1:nrow(gsea.result)) {
  gs.name <- gsea.result$NAME[i]
  gs.genes.file <- paste0(opt$'gene-set-dir', "/", gs.name, ".xls")
  if (!file.exists(gs.genes.file)) {
    #print(paste("Output file for gene set ", gs.name, "not found. SKIPPED!"))
    skipped <- skipped + 1
    next
  }
  #print(paste("Adding ", i, "of", nrow(gsea.result), "gene sets:", gs.name))
  gs.genes <- read.delim(gs.genes.file, row.names = 1, stringsAsFactors = FALSE)
  gs.genes <- gs.genes[gs.genes$CORE.ENRICHMENT == "Yes",]
  
  if (opt$phenotype == 'up') {
    gs.genes$rank <- gs.genes$RANK.IN.GENE.LIST
  } else {
    gs.genes$rank <- opt$'num-genes' - gs.genes$RANK.IN.GENE.LIST + 1
  }
  gs.genes <- gs.genes[order(gs.genes$rank),]
  gs.genes$CORE_GENES <- paste0(gs.genes$PROBE, "(", gs.genes$rank, ")")
  collapsed <- paste(gs.genes$CORE_GENES, collapse=", ")
  out$CORE_GENES[out$NAME==gs.name] <- collapsed

  processed <- processed + 1
}

print(paste("Processed", processed, "gene sets."))
if (skipped > 0) print(paste("WARNING: Skipped", skipped, "gene sets because output file was not found. Consider increasing '-plot_top_x' value when running GSEA."))

# sort
out <- out[order(abs(out$NES), decreasing = TRUE),]

# write output
names(out)[2] <- "DESCRIPTION"
write.table(out[,c("NAME", "DESCRIPTION", "CATEGORY", "SIZE", "ES", "NES", "NOM p-val", "FDR q-val", "FWER p-val", "RANK AT MAX", "LEADING EDGE", "CORE_GENES")], opt$'output-file', quote=FALSE, sep="\t",  na="", row.names=FALSE, col.names=TRUE)

# ---
# experimental: determine overlaps of core genes between gene sets to reduce redundancy in output
# however, this did not work well because overlap is not that large for most gene sets (Jaccard index < 0.5) :-(

#name.sets <- unique(combined$geneset)
#num.sets <- length(name.sets)
#m <- matrix(as.numeric(0), nrow=num.sets, ncol=num.sets, dimnames=list(name.sets, name.sets))
#genes.sets <- lapply(1:num.sets, function(x) { combined$gene[combined$geneset==name.sets[x]] })
#for (i in 1:num.sets) {
#  for (j in 1:num.sets) {
#    m[i, j] = length(intersect(genes.sets[[i]], genes.sets[[j]])) / length(union(genes.sets[[i]], genes.sets[[j]]))
#  }
#}

#library(gplots)
#heatmap.2(m, trace="none", scale="none")
#m[m>0.8 & m < 1]
#which(m == 0.75, arr.ind = TRUE)

