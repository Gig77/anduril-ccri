library(componentSkeleton)

execute <- function(cf) {

	# debug
	#rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/fikret/results/anduril/execute/degBoxplotUp_DTCvsMNC-geneBoxplot/_command")

	instance.name <- get.metadata(cf, 'instanceName')
	
	# ----
	
	expr <- Matrix.read(get.input(cf, 'expr'))
	sampleGroups <- CSV.read(get.input(cf, "sampleGroups"))
	annotation <- CSV.read(get.input(cf, "annotation"))
	
	# read sample groups
	sample2group <- data.frame(sample=character(0), group=character(0))
	for(g in sampleGroups$ID) {
	  members <- unlist(strsplit(sampleGroups[sampleGroups[,'ID']==g,'Members'],','))
	  sample2group <- rbind(sample2group, data.frame(sample=members, group=g, stringsAsFactors = F))
	}
	
	# only included groups
	includeGroups <- get.parameter(cf, 'includeGroups', 'string')
	if (includeGroups != "") {
	  sample2group <- sample2group[sample2group$group %in% unlist(strsplit(includeGroups, ',')),]
	  if (nrow(sample2group) == 0) {
		  stop(sprintf("ERROR: No sample found belonging to one of specified groups: %s", includeGroups))
	  }
	}
	
	# translate ensembl gene ids to HGNC symbols
	names(annotation)[names(annotation)=="Ensembl Gene ID"] <- "Ensembl"
	names(annotation)[names(annotation)=="Associated Gene Name"] <- "HGNC"
	
	hgncs <- get.parameter(cf, 'hgnc', 'string')
	
	genesViaParameter <- !is.null(hgncs) && !is.na(hgncs) && hgncs != ""
	if (genesViaParameter) {  # gene list (HGNC symbols) provided via parameter?
	  genes <- data.frame(HGNC=unlist(strsplit(hgncs, ',')), stringsAsFactors = F)
	  genes <- merge(genes, annotation[,c("Ensembl", "HGNC")], all.x=T)
	  genes$Ensembl[is.na(genes$Ensembl)] <- genes$HGNC[is.na(genes$Ensembl)]
	} else {  # no: use genes (ensembl ids) provided via table input port
	  geneIds <- CSV.read(get.input(cf, "geneIds"))
	  names(geneIds)[1] <- "Ensembl"
	  genes <- merge(subset(geneIds, select="Ensembl"), annotation[,c("Ensembl", "HGNC")], all.x=T)
	  name.missing <- is.na(genes$HGNC) | genes$HGNC == ""
	  genes$HGNC[name.missing] <- genes$Ensembl[name.missing]
	}
	genes <- merge(genes, expr, by.x="Ensembl", by.y="row.names")

	# prepare document
	tex <- character(0)
	if (nrow(genes) > 0) tex <- c(tex, '\\clearpage')
	
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
	
	if (nrow(genes) == 0) {
	  tex <- c(tex, "No significant DEGs found.")
	} else {
	  
  	# melt dataframe into shape suitable for plotting and add sample groups
  	library(reshape)
  	gexpr <- melt(genes, id.vars=c("Ensembl", "HGNC"))
  	names(gexpr)[names(gexpr)=="variable"] <- "sample"
  	gexpr$sample <- as.character(gexpr$samp)
  	gexpr <- merge(gexpr, sample2group, all.X=T)
  	
  	# before plotting, preserve input gene order in lattice output (note: factor levels determines panel order)
  	# original order might have been destroyed by merge
  	if (genesViaParameter) {
  	  genes <- genes[match(unlist(strsplit(hgncs, ',')), genes$HGNC),] 
  	} else {
  	  matches <- match(geneIds$Ensembl, genes$Ensembl)
  	  matches <- matches[!is.na(matches)]
  	  genes <- genes[matches,] # original order was destroyed by merge, so let's restore it first
  	}                        	
  	gexpr$HGNC <- factor(as.character(gexpr$HGNC), levels=genes$HGNC)
  
  	if (includeGroups != "") {
  		gexpr$group <- factor(as.character(gexpr$group), levels=unlist(strsplit(includeGroups, ',')))
  		if(sum(is.na(gexpr$group)) > 0) stop(sprintf("Incomplete assignment of samples to groups due to invalid or missing sample group(s) in parameter 'includeGroups': %s", includeGroups))
  	}
  	
  	# plot page by page
  	library(lattice)
  	
  	nRow <- get.parameter(cf, 'nRow', 'int')
  	nCol <- get.parameter(cf, 'nCol', 'int')
  		
  	splits <- split(levels(gexpr$HGNC), (0:(length(levels(gexpr$HGNC))-1)) %/% (nRow*nCol) + 1)
  
  	width <- get.parameter(cf, 'width', 'float')
  	height <- get.parameter(cf, 'height', 'float')
  	label.outliers <- get.parameter(cf, 'labelOutliers', 'boolean')
  	cex.dot <- get.parameter(cf, 'cexDot', 'float')
  	cex.sample.label <- get.parameter(cf, 'cexSampleLabel', 'float')
	  cex.group.label <- get.parameter(cf, 'cexGroupLabel', 'float')
	
  	rowheight <- (height-2)/nRow
  	
  	for(pageno in 1:length(splits)) {
  	  gexpr.thispage <- gexpr[gexpr$HGNC %in% splits[[pageno]],]
  	  nRow.thispage <- min(nRow, ((length(unique(gexpr.thispage$HGNC))-1) %/% nCol) + 1)
  	  
  	  # get output file name
  	  plot.file <- sprintf('geneboxplot-%s-figure%d.pdf', instance.name, pageno)
  	  
  	  pdf(file.path(out.dir, plot.file), height=rowheight*nRow.thispage+2, width=width)
  	  print(bwplot(value~group | HGNC, data=gexpr.thispage, 
  	               layout=c(nCol,nRow.thispage),
  	               par.strip.text=list(cex=0.5),
  	               notch=FALSE,
  	               as.table=TRUE,
  	               ylab="DESeq2 normalized expression",
  	               scales=list(x=list(rot=90,cex=cex.group.label)),
  	               par.settings = list(box.umbrella=list(col="black"), box.rectangle = list(col="black")), 
  	               panel=function(x,y,...,subscripts){
  	                 panel.grid()
  	                 bw <- panel.bwplot(x,y,pch="|",do.out=FALSE, ...)
  	                 panel.stripplot(x,y,jitter.data=TRUE,factor=0.8,pch=19,cex=cex.dot, ...)
  	                 if (label.outliers) {
  	                   whisker.up <- tapply(y, factor(x), FUN=function(d) { boxplot.stats(d)$stats[5]})
  	                   whisker.dn <- tapply(y, factor(x), FUN=function(d) { boxplot.stats(d)$stats[1]})
  	                   lab <- as.character(gexpr.thispage$sample[subscripts])
  	                   lab[y <= whisker.up[x] & y >= whisker.dn[x]] <- ""
  	                   panel.text(x, y, labels=lab, cex=cex.sample.label)
  	                 }
  	               }))
  	  dev.off()
  	  
  	  # generate latex string
  	  caption <- get.parameter(cf, 'caption', 'string')
  	  if (length(splits) > 1) {
  	    caption <- paste0(caption, " Part ", pageno, " of ", length(splits), ".")
  	  }
  	  tex <- c(tex, latex.figure(plot.file, caption=caption))
  	}
  }

	latex.write.main(cf, 'document', tex)
	
	return(0)
}

main(execute)

