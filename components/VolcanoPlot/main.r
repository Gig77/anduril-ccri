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
	caption <- get.parameter(cf, 'caption', 'string')
	sampleSize <- get.parameter(cf, 'sampleSize', 'int')
	cexLabel <- get.parameter(cf, 'cexLabel', 'float')
	
	# prepare data
	res <- data.frame(id=expr[,1], log2FoldChange=expr[,fccol], pvalue=pmax(expr[,pcol], 1e-200), padj=expr[,qcol])

	if (get.input(cf, "geneNames") != "") {
	  geneNames <- CSV.read(get.input(cf, "geneNames"))
	  res <- merge(res, geneNames[,c(1,geneNameCol)], by.x=1, by.y=1, all.x=T)
	  res <- res[!duplicated(res$id),]
	  names(res)[5] <- "Name"
	} else {
	  res$Name <- res[,1]
	}

	# subsample uninteresting fraction of genes to reduce plot size
	res.interesting   <- subset(res, padj <= sigthresh | abs(log2FoldChange) >= lfcthresh)
	res.uninteresting <- subset(res, padj > sigthresh & abs(log2FoldChange) < lfcthresh)
	res.uninteresting.sample <- res.uninteresting[sample(1:nrow(res.uninteresting), min(sampleSize, nrow(res.uninteresting)), replace = F),]
  	res.sample <- rbind(res.interesting, res.uninteresting.sample)
  
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
	with(subset(res.sample, padj >  sigthresh & abs(log2FoldChange) <  lfcthresh), plot(log2FoldChange, -log10(pvalue), pch=20, col="lightgray", xlim=range(res.sample$log2FoldChange), ylim=range(-log10(res.sample$pvalue)), main=""))
  with(subset(res.sample, padj <= sigthresh & abs(log2FoldChange) <  lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="black"))
	with(subset(res.sample, padj >  sigthresh & abs(log2FoldChange) >= lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
	with(subset(res.sample, padj <= sigthresh & abs(log2FoldChange) >= lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  if (labelTopN > 0) {
    library(calibrate)
	res.sample.sorted <- res.sample[order(res.sample$padj, -abs(res.sample$log2FoldChange)),]
	res.sample.sorted <- res.sample.sorted[!is.na(res.sample.sorted$padj),]
	with(res.sample.sorted[1:min(labelTopN, nrow(res.sample.sorted)),], textxy(log2FoldChange, -log10(pvalue), labs=Name, cex=cexLabel, offset=0.7))
	caption <- paste0(caption, " The top-", labelTopN, " most significant genes are labeled.")
}
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<=",sigthresh,sep=""), paste("|LogFC|>=",lfcthresh,sep=""), "both"), pch=20, col=c("black","blue","red"))
  dev.off()
	  
  # generate latex string
  tex <- c(tex, latex.figure(plot.file, caption=caption))
	latex.write.main(cf, 'document', tex)

	return(0)
}

main(execute)
