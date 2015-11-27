library(componentSkeleton)

execute <- function(cf) {

	# enable for debugging
	#cf <- parse.command.file("/mnt/projects/martin/results/anduril/execute/qcReport-samplesClusterHeatmap/_command")
	
	# Meta data
	instance.name <- get.metadata(cf, 'instanceName')	

	# Inputs
	countMatrix  <- Matrix.read(get.input(cf, "countMatrix"))
	samples      <- CSV.read(get.input(cf, 'samples'))
	
	# Params
	caption <- get.parameter(cf, 'caption',    type = 'string')
	cexRow <- get.parameter(cf, 'cexRow',    type = 'float')
	section.title <- get.parameter(cf, 'sectionTitle', type = 'string')
	section.type <- get.parameter(cf, 'sectionType', type = 'string')
	do.rlog <- get.parameter(cf, 'rlog', type = 'boolean')
	do.vst <- get.parameter(cf, 'vst', type = 'boolean')
	do.voom <- get.parameter(cf, 'voom', type = 'boolean')
			
	library("DESeq2")
	library("RColorBrewer")
	library("gplots")
	
	rownames(samples) <- samples$Alias
	samples <- samples[samples$Alias %in% colnames(countMatrix),]
	countMatrix <- countMatrix[,rownames(samples)]
	cds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = samples, design = ~1)	
	expressed <- rowSums(counts(cds)) >= 10
	caption <- paste(caption, "The clustering is based on", sum(expressed), "expressed genes.", sum(!expressed), "genes were excluded because they had less than ten reads in total across all samples.")
	hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

	report.dir <- get.output(cf, 'report')
	dir.create(report.dir, recursive=TRUE)
	
	# rlog
	if (do.rlog) {
		rld <- rlog(cds[expressed,])
		CSV.write(get.output(cf, 'rlog'), assay(rld), first.cell = "Ensembl")
		#dists.rld <- as.matrix(dist(t(assay(rld))))
		dists.rld <- as.matrix(as.dist(1-cor(assay(rld), method="spearman"))) # spearman correlation
		
		plot.file.rld <- sprintf('%s-heatmap-rlog.pdf', instance.name)
		pdf(file.path(report.dir, plot.file.rld))
		heatmap.2(dists.rld, trace="none", scale="none", col=rev(hmcol), margin=c(7, 7), cexRow=cexRow, cexCol=cexRow, key.title="")
		dev.off()
	} else {
	  file.create(get.output(cf, 'rlog'))
	}
	
	# VST
	if (do.vst) {
		vst <- varianceStabilizingTransformation(cds[expressed,])
		CSV.write(get.output(cf, 'vst'), assay(vst), first.cell = "Ensembl")
		dists.vst <- as.matrix(as.dist(1-cor(assay(vst), method="spearman"))) # spearman correlation
		
		plot.file.vst <- sprintf('%s-heatmap-vst.pdf', instance.name)
		pdf(file.path(report.dir, plot.file.vst))
		heatmap.2(dists.vst, trace="none", scale="none", col=rev(hmcol), margin=c(7, 7), cexRow=cexRow, cexCol=cexRow, key.title="")
		dev.off()
	} else {
	  file.create(get.output(cf, 'vst'))
	}
	
	# voom
	if (do.voom) {
		library(edgeR)
		dge <- DGEList(counts=counts(cds[expressed,]))
		dge.norm <- calcNormFactors(dge, method="TMM")
		y <- voom(dge.norm)
		#dists.voom <- as.matrix(dist(t(y$E)))
		CSV.write(get.output(cf, 'voom'), y$E, first.cell = "Ensembl")
		dists.voom <- as.matrix(as.dist(1-cor(y$E, method="spearman"))) # spearman correlation
		
		plot.file.voom <- sprintf('%s-heatmap-voom.pdf', instance.name)
		pdf(file.path(report.dir, plot.file.voom))
		heatmap.2(dists.voom, trace="none", scale="none", col=rev(hmcol), margin=c(7, 7), cexRow=cexRow, cexCol=cexRow, key.title="")
		dev.off()
	} else {
	  file.create(get.output(cf, 'voom'))
	}
	
	# prepare Latex document
	tex <- character()
	tex <- c(tex, '\\clearpage')
	if (nchar(section.type) > 0) {
		section.type=paste('\\',section.type,sep='') # if section type is defined, add escape in the beginning
	}
	if (nchar(section.title) > 0) {
		tex <- c(tex, sprintf('%s{%s}\\label{%s}', section.type, section.title,	instance.name))
	}
	
	if (do.rlog) {
		tex <- c(tex, latex.figure(plot.file.rld, caption=paste0(caption, " Before clustering, read counts were transformed by DESeq2 regularized log (rlog) transformation and sample distances were computed using the Spearman correlation coefficient.")))		    			    
	}
	if (do.vst) {	
		tex <- c(tex, latex.figure(plot.file.vst, caption=paste0(caption, " Before clustering, read counts were transformed by DESeq2 variance stabilizing transformation (VST) and sample distances were computed using the Spearman correlation coefficient.")))		    			    
	}
	if (do.voom) {
		tex <- c(tex, latex.figure(plot.file.voom, caption=paste0(caption, " Before clustering, read counts were transformed by edgeR Voom and sample distances were computed using the Spearman correlation coefficient.")))		    			    
	}
	latex.write.main(cf, 'report', tex)
	
	return(0)
}

main(execute)

