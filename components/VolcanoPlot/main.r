library(componentSkeleton)

execute <- function(cf) {

	# debug
	#rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/helena_veronika/results/anduril/execute/degVolcanoPlot_oeERvsRHD/_command")

	instance.name <- get.metadata(cf, 'instanceName')
	
	# get input data and parameters
	expr <- CSV.read(get.input(cf, 'expr'))

	sigthresh <- get.parameter(cf, 'sigthresh', 'float')
	lfcthresh <- get.parameter(cf, 'lfcthresh', 'float')
	fccol <- get.parameter(cf, 'fccol', 'int')
	pcol <- get.parameter(cf, 'pcol', 'int')
	qcol <- get.parameter(cf, 'qcol', 'int')
	geneNameCol <- get.parameter(cf, 'geneNameCol', 'int')
	labelTopN <- get.parameter(cf, 'labelTopN', 'int')
	legendpos <- get.parameter(cf, 'legendpos', 'string')
	excludeGenes <- get.parameter(cf, 'excludeGenes', 'string')
	caption <- get.parameter(cf, 'caption', 'string')
	sampleSize <- get.parameter(cf, 'sampleSize', 'int')
	cexLabel <- get.parameter(cf, 'cexLabel', 'float')
	cexPoint <- get.parameter(cf, 'cexPoint', 'float')
	minP <- get.parameter(cf, 'minP', 'float')
	
	# prepare data
	res <- data.frame(id=expr[,1], log2FoldChange=expr[,fccol], pvalue=pmax(expr[,pcol], minP), padj=expr[,qcol], stringsAsFactors = F)
	res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange),]

	if (get.input(cf, "geneNames") != "") {
	  geneNames <- CSV.read(get.input(cf, "geneNames"))
	  res <- merge(res, geneNames[,c(1,geneNameCol)], by.x=1, by.y=1, all.x=T)
	  res <- res[!duplicated(res$id),]
	  names(res)[5] <- "Name"
	  res$Name[is.na(res$Name) | res$Name == ""] <- res$id[is.na(res$Name) | res$Name == ""]
	} else {
	  res$Name <- res[,1]
	}

	# subsample uninteresting fraction of genes to reduce plot size
	res.interesting   <- subset(res, padj <= sigthresh | abs(log2FoldChange) >= lfcthresh)
	res.uninteresting <- subset(res, padj > sigthresh & abs(log2FoldChange) < lfcthresh)
	res.uninteresting.sample <- res.uninteresting[sample(1:nrow(res.uninteresting), min(sampleSize, nrow(res.uninteresting)), replace = F),]
  	res.sample <- rbind(res.interesting, res.uninteresting.sample)
	
	if (excludeGenes != "") {
		exclude <- unlist(strsplit(excludeGenes,','))
		res.sample <- res.sample[!res.sample$Name %in% exclude,]
	}
  
	# prepare document
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

	plot.file <- sprintf('volcanoplot-%s.pdf', instance.name)
	
	# plot
	pdf(file.path(out.dir, plot.file))
	with(subset(res.sample, padj >  sigthresh & abs(log2FoldChange) <  lfcthresh), plot(log2FoldChange, -log10(pvalue), pch=20, col="lightgray", xlim=c(min(res.sample$log2FoldChange)-1,max(res.sample$log2FoldChange)+1), ylim=range(-log10(res.sample$pvalue)), cex=cexPoint, main=""))
    with(subset(res.sample, padj <= sigthresh & abs(log2FoldChange) <  lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="black", cex=cexPoint))
	with(subset(res.sample, padj >  sigthresh & abs(log2FoldChange) >= lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="blue", cex=cexPoint))
	with(subset(res.sample, padj <= sigthresh & abs(log2FoldChange) >= lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="red", cex=cexPoint))

	# label points
	if (labelTopN > 0) {
    library(calibrate)
    
  	res.sample.sorted <- res.sample[order(res.sample$padj),]
  	res.sample.sorted <- res.sample.sorted[!is.na(res.sample.sorted$padj) & res.sample.sorted$log2FoldChange > 0,]
	label.up.sig <- res.sample.sorted$id[1:min(labelTopN, nrow(res.sample.sorted))]

	res.sample.sorted <- res.sample[order(-res.sample$log2FoldChange),]
	res.sample.sorted <- res.sample.sorted[!is.na(res.sample.sorted$padj) & res.sample.sorted$log2FoldChange > 0,]
	label.up.fc <- res.sample.sorted$id[1:min(labelTopN, nrow(res.sample.sorted))]
	
	with(res.sample.sorted[res.sample.sorted$id %in% c(label.up.sig, label.up.fc),], textxy(log2FoldChange, -log10(pvalue), labs=Name, cex=cexLabel, offset=0.7))

  	res.sample.sorted <- res.sample[order(res.sample$padj),]
  	res.sample.sorted <- res.sample.sorted[!is.na(res.sample.sorted$padj) & res.sample.sorted$log2FoldChange < 0,]
	label.dn.sig <- res.sample.sorted$id[1:min(labelTopN, nrow(res.sample.sorted))]

	res.sample.sorted <- res.sample[order(res.sample$log2FoldChange),]
	res.sample.sorted <- res.sample.sorted[!is.na(res.sample.sorted$padj) & res.sample.sorted$log2FoldChange < 0,]
	label.dn.fc <- res.sample.sorted$id[1:min(labelTopN, nrow(res.sample.sorted))]
	
  	with(res.sample.sorted[res.sample.sorted$id %in% c(label.dn.sig, label.dn.fc),], textxy(log2FoldChange, -log10(pvalue), labs=Name, cex=cexLabel, offset=0.7))
  	
  	caption <- paste0(caption, " The top-", labelTopN, " most significant up- and down-regulated genes (either in terms of q-value or fold-change) are labeled. P-values were capped at ", sprintf("%.1g", minP), ".")
	if (excludeGenes != "") {
		caption <- paste0(caption, " The following genes are not shown: ", excludeGenes, ".")
	}
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<=",sigthresh,sep=""), paste("|LogFC|>=",lfcthresh,sep=""), "both"), pch=20, col=c("black","blue","red"))
  dev.off()
	  
  # generate latex string
  tex <- c(tex, latex.figure(plot.file, caption=caption))
	latex.write.main(cf, 'document', tex)

	return(0)
}

main(execute)
