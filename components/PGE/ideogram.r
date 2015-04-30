library(ggbio)
library(GenomicRanges)
library(optparse)
library(dplyr)

option_list <- list(
		make_option("--query", type="character", help="Query genes used for enrichement"),
		make_option("--reference", type="character", help="List of all genes in genome"),
		make_option("--regions", type="character", help="Enriched regions"),
		make_option("--maxRegionsPerChr", type="integer", help="Enriched regions"),
		make_option("--output", type="character", help="Name of PDF output file")
)
opt <- parse_args(OptionParser(option_list=option_list))
#opt <- data.frame(query="/mnt/synology/data/christian/iamp/results/current/anduril/execute/degs_ERvsNonER/table.csv", reference="/mnt/synology/data/christian/iamp/scripts/pge/pge_release_1.01/ensembl75.txt", regions="/mnt/synology/data/christian/iamp/results/current/anduril/execute/pge_ERvsNonER-pge/enrichedRegions.csv", maxRegionsPerChr=20, output="/mnt/synology/data/christian/iamp/results/current/ideogram_test.pdf", stringsAsFactors=F, check.names=F)

reference <- read.delim(opt$reference, header=F, stringsAsFactors=F)
colnames(reference) <- c("ids", "chr", "start", "end")

query <- read.delim(opt$query, stringsAsFactors=F)
query.gr <- merge(query, reference)
query.gr <- GRanges(seqnames=paste0("chr", query.gr$chr), ranges=IRanges(start=query.gr$start, end=query.gr$end))

regions <- read.delim(opt$regions, stringsAsFactors=F)

if (nrow(regions) > 0) {
	regions <- group_by(regions, chr) %>% top_n(opt$maxRegionsPerChr, desc(pvalue)) # keep only top n regions per chromosome
	regions <- regions[order(regions$end-regions$start, decreasing=T),] # plot largest regions closest to chromosome
	regions.gr <- GRanges(seqnames=paste0("chr", regions$chr), ranges=IRanges(start=pmin(regions$start, regions$end), end=pmax(regions$start, regions$end)), pvalue=regions$pvalue)
	
	# alpha value (transparency) propotional to p-value, scale between 0 and 1, assign 1 (least transparency) for smallest p-value
	regions.gr$alpha <- pmin(1, (-log(regions.gr$pvalue)-min(-log(regions.gr$pvalue))) / (max(-log(regions.gr$pvalue))-min(-log(regions.gr$pvalue)))+0.2)
				
	# stack algorithm; place non-overlapping regions next to each other
	o <- findOverlaps(regions.gr, regions.gr)
	o <- o[o@queryHits != o@subjectHits] # remove self-hits
	regions.gr$level <- as.integer(0)
	placed <- list()
	regions.gr$level[1] <- 1
	placed[[1]] <- 1
	if (length(regions.gr) > 1) {
		for (i in 2:length(regions.gr)) {
			
			level <- 1
			while (level <= length(placed) && length(intersect(o@subjectHits[o@queryHits==i], placed[[level]])) > 0) { 
				level <- level + 1 
			}
			regions.gr$level[i] <- level
			if (level > length(placed)) {
				placed[[level]] <- i
			} else {
				placed[[level]] <- c(placed[[level]], i)
			}
		}
	}
}

# draw ideogram
data(hg19IdeogramCyto, package = "biovizBase")
hg19 <- keepSeqlevels(hg19IdeogramCyto, paste0("chr", c(1:22, "X", "Y")))

pdf(opt$output)
p <- ggplot(hg19) + theme(legend.position = "none", axis.title.x=element_blank(), axis.text.x=element_blank()) 
p <- p + layout_karyogram(cytoband = TRUE) 
p <- p + layout_karyogram(query.gr, geom = "rect", aes(xmin = start, xmax = end, ymin = 12, ymax = 20), fill = "red", size=0)
if (nrow(regions) > 0) {
	p <- p + layout_karyogram(regions.gr, geom = "rect", aes(xmin = start, xmax = end, ymin = -((level-1)*2+2), ymax = -((level-1)*2+4), alpha=alpha), fill = "blue", size=0)
}
print(p)
dev.off()

