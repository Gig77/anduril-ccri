library(componentSkeleton)

execute <- function(cf) {

	# debug
	#rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/iamp/results/anduril/execute/qcReport-libraryBatches/_command")
  #stop("HERE!")
  
  instance.name <- get.metadata(cf, 'instanceName')
  
	# read count metrics
  readStats <- CSV.read(get.input(cf, "readStats"))
  readStats$pct_rrna <- (readStats$exonic-readStats$'non-rRNA') / readStats$exonic
	readStats$pct_rna <- (readStats$'non-rRNA'-readStats$protein) / readStats$'non-rRNA'
	readStats$pct_intronic <- (readStats$mapped-readStats$exonic) / readStats$mapped
	readStats$pct_dup <- (readStats$'uniquely mapped'-readStats$'non-duplicates') / readStats$'uniquely mapped'
	
	# 5' 3' coverage bias
	cov <- as.data.frame(Matrix.read(get.input(cf, 'geneBodyCoverages')))
	cov$five_three_prime_cov_bias <- rowSums(cov[,1:25])/rowSums(cov[,75:100])

	# GC bias
	#---
	# NOTE: GC expression ratio of housekeeping genes is currently not used for sample clustering
	# because we might introduce and unknown amount of true biological variation due to true differential
	# expression of these genes; if included in PCA, we get an almost perfect separation of samples
	# by biological subtype, which is suspicious...
	#---
	#gcRatio <- getGCRatio()
	
  # collect all metrics in dataframe	
	batch <- merge(readStats[,c("sample", "pct_rrna", "pct_rna", "pct_intronic", "pct_dup")], subset(cov, select=c("five_three_prime_cov_bias")), by.x='sample', by.y="row.names")
	#batch <- merge(batch, subset(gcRatio, select=c("gcRatio")), by.x='sample', by.y="row.names")
	
	# perform PCA on ratios
	pca <- prcomp(batch[,c("pct_rrna", "pct_intronic", "pct_dup", "five_three_prime_cov_bias")], center = TRUE, scale. = TRUE) 
	
	# transform PCs into categorial values
	batch$PC1 <- pca$x[,1]
	batch$PC1_categorial <- NA
	batch$PC1_categorial[pca$x[,1] <= quantile(pca$x[,1], 0.3)] <- "low"
	batch$PC1_categorial[pca$x[,1] >= quantile(pca$x[,1], 0.7)] <- "high"
	batch$PC1_categorial[is.na(batch$PC1_categorial)] <- "norm"
	
	batch$PC2 <- pca$x[,2]
	batch$PC2_categorial <- NA
	batch$PC2_categorial[pca$x[,2] <= quantile(pca$x[,2], 0.3)] <- "low"
	batch$PC2_categorial[pca$x[,2] >= quantile(pca$x[,2], 0.7)] <- "high"
	batch$PC2_categorial[is.na(batch$PC2_categorial)] <- "norm"
	
	# add sample groups
	batch$group <- NA
	groups <- componentSkeleton::SampleGroupTable.read(get.input(cf, "groups"))
	for (gid in componentSkeleton::SampleGroupTable.get.groups(groups)) {
	  for(s in componentSkeleton::SampleGroupTable.get.source.groups(groups, gid)) {
	    batch$group[batch$sample==s] <- gid
	  }
	}
	batch$group <- factor(batch$group, levels=groups$ID)
	
	CSV.write(get.output(cf, 'batches'), batch)
	
	# PCA plot
	
	tex <- character(0)
	tex <- c(tex, '\\clearpage')
	
	out.dir <- get.output(cf, 'document')
	dir.create(out.dir, recursive=TRUE)
	
	section.title <- get.parameter(cf, 'sectionTitle')
	section.type <- get.parameter(cf, 'sectionType')
	if (nchar(section.type) > 0) {
	  section.type=paste('\\',section.type,sep='') # if section type is defined, add escape in the beginning
	}
	if (nchar(section.title) > 0) {
	  tex <- c(tex, sprintf('%s{%s}\\label{%s}', section.type, section.title,	instance.name))
	}
	
	plot.file <- sprintf('pcaplot-%s.pdf', instance.name)
	pdf(file.path(out.dir, plot.file))
	
	par(mar=c(4,4,4,7))
	par(xpd=TRUE)
	plot(pca$x[,1], pca$x[,2], col="white", xlab="PC1", ylab="PC2")
	text(pca$x[,1], pca$x[,2], batch$sample, col=as.numeric(batch$group), cex=0.8)
	legend(x=par("usr")[2], y=max(par("usr")[3:4]), levels(batch$group), col=1:length(levels(batch$group)), pch=1)

	dev.off()
	
	# generate latex string
	caption <- get.parameter(cf, 'caption', 'string')
	tex <- c(tex, latex.figure(plot.file, caption=caption))

	latex.write.main(cf, 'document', tex)

	return(0)
}

getGCRatio <- function() {
  # GC bias: get GC content and length of (longest) transcript per gene
  
  # get genes expressed in our experiment
  expr <- Matrix.read(get.input(cf, "expr"))
  exprGenes <- rownames(expr)[rowSums(expr) >= 10]
  
  # look only at housekeeping genes b/c they should be more stably expressed across conditions
  # library(BiocInstaller) ; biocLite("tweeDEseqCountData")  # 40 Mb
  library(tweeDEseqCountData)
  data(hkGenes)
  exprHkGenes <- intersect(exprGenes, hkGenes)
  
  library("biomaRt")
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
  result <- select(mart, keys=exprHkGenes, keytype="ensembl_gene_id", columns=c("ensembl_gene_id", "percentage_gc_content", "transcript_length"))
  result.longest <- result[ave(result$transcript_length, result$ensembl_gene_id, FUN=max) == result$transcript_length,]
  feature <- data.frame(row.names=result.longest$ensembl_gene_id, gc=result.longest$percentage_gc_content, length=result.longest$transcript_length)
  
  gcpoor <- rownames(feature)[order(feature$gc)][1:100]
  gcrich <- rownames(feature)[order(feature$gc, decreasing=TRUE)][1:100]
  
  gcRatio <- as.data.frame(as.matrix(colSums(expr[gcpoor,])/colSums(expr[gcrich,])))
  names(gcRatio) <- "gcRatio"
  
  retrun(gcRatio)
}

main(execute)

