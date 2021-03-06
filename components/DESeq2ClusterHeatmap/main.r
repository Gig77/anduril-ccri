library(componentSkeleton)

execute <- function(cf) {

	# enable for debugging
	# rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/fikret/results/anduril/execute/qcReport-samplesClusterHeatmap/_command")
	
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
	annotations <- get.parameter(cf, 'annotations', type = 'string')
	legendYOffset <- get.parameter(cf, 'legendYOffset', type = 'float')
	legendXOffset <- get.parameter(cf, 'legendXOffset', type = 'float')
	margins <- get.parameter(cf, 'margins', type = 'string')
	method <- get.parameter(cf, 'method', type = 'string')
	minSumReads <- get.parameter(cf, 'minSumReads', type = 'int')
	
	library("DESeq2")
	library("RColorBrewer")
	library("gplots")
	library("heatmap3")
	
	rownames(samples) <- samples$Alias
	samples <- samples[samples$Alias %in% colnames(countMatrix),]
	countMatrix <- countMatrix[,rownames(samples)]
	cds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = samples, design = ~1)	
	expressed <- rowSums(counts(cds)) >= minSumReads
	caption <- paste(caption, "The clustering is based on", sum(expressed), "expressed genes.", sum(!expressed), "genes were excluded because they had less than", minSumReads, " reads in total across all samples.")

	report.dir <- get.output(cf, 'report')
	dir.create(report.dir, recursive=TRUE)
	
	plotHeatmap <- function(matr, samples, annotations, cexRow) {
	  hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
	  noLegend <- function() showLegend(legend=c(""), col="white")
	  
	  # setup annotation colors drawn on top of heatmap
	  if (annotations != "") {
	    annotations <- unlist(strsplit(annotations, ","))
	    ann.factors <- list()
	    ann.palettes <- list()
	    ann.colors <- NULL
	    for(a in annotations) {
	      ann.cur <- as.character(samples[,a])
	      ann.cur[is.na(ann.cur) | ann.cur == ""] <- "n/a"
	      levels <- unique(ann.cur)[unique(ann.cur) != "n/a"]
	      if(sum(ann.cur == "n/a") > 0) levels <- c(levels, "n/a")
	      v <- factor(ann.cur, levels=levels)
	      
	      # choose set of distinguishable colors for sample annotations
	      # see also http://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
	      if (length(levels) <= 11) {
	        ann.palettes[[a]] <- sample(c("#FFFF99", "#386CB0", "#BF5B17", "#666666", "#7FC97F", "#BEAED4", "#FF7F00", "#FFFF33", "#999999", "#F2F2F2", "#000000"), length(levels))
	      } else {
	        set.seed(4711)
	        ann.palettes[[a]] <- sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert=T)], length(levels))
	      }
	      if (is.null(ann.colors)) {
	        ann.colors <- cbind(ann.palettes[[a]][v])
	      } else {
	        ann.colors <- cbind(ann.colors, ann.palettes[[a]][v]) 
	      }
	      ann.factors[[a]] <- v
	    } 
	    dimnames(ann.colors)[[2]] <- annotations
	    
	    # draw heatmap including annotations
	    heatmap3(matr, col=rev(hmcol), scale="none", cexRow=cexRow, cexCol=cexRow, method=method, margins=as.numeric(unlist(strsplit(margins, ","))), ColSideColors = ann.colors, legendfun = noLegend)
	    
	    # draw legends into right margin
	    twidth <- max(strwidth(names(ann.factors))) # to have equally sized legends and thus left-aligned boxes
	    y = 1-legendYOffset
	    for (a in names(ann.factors)) {
	      v <- levels(ann.factors[[a]])
	      x <- legend(1-legendXOffset, y, v, fill=ann.palettes[[a]], bty="n", cex=0.5, title=a, xjust = 0, title.adj = 0, text.width=twidth, y.intersp = 0.8, xpd=TRUE)
	      y = y - 0.027 * length(v) - 0.02
	    }
	  } else {
	    # draw heatmap without annotations on top
	    heatmap3(matr, col=rev(hmcol), scale="none", cexRow=cexRow, cexCol=cexRow, method=method, margins=as.numeric(unlist(strsplit(margins, ","))), legendfun = noLegend)
	  }
	}
	
	# rlog
	if (do.rlog) {
		rld <- rlog(cds[expressed,])
		CSV.write(get.output(cf, 'rlog'), assay(rld), first.cell = "Ensembl")
		#dists.rld <- as.matrix(dist(t(assay(rld))))
		dists.rld <- as.matrix(as.dist(1-cor(assay(rld), method="spearman"))) # spearman correlation
		
		plot.file.rld <- sprintf('%s-heatmap-rlog.pdf', instance.name)
		pdf(file.path(report.dir, plot.file.rld))
		plotHeatmap(matr=dists.rld, samples=samples, annotations=annotations, cexRow=cexRow)
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
		plotHeatmap(matr=dists.vst, samples=samples, annotations=annotations, cexRow=cexRow)
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
		plotHeatmap(matr=dists.voom, samples=samples, annotations=annotations, cexRow=cexRow)
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

