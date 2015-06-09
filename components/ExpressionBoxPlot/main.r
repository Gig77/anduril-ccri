library(componentSkeleton)

execute <- function(cf) {

	# debug
	#rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/iamp/results/anduril/execute/degBoxplotUp_ERvsPreB-geneBoxplot/_command")

	instance.name <- get.metadata(cf, 'instanceName')
	
	# ----
	
	expr <- Matrix.read(get.input(cf, 'expr'))
	groups <- CSV.read(get.input(cf, "groups"))
	annotation <- CSV.read(get.input(cf, "annotation"))
	
	# read sample groups
	sample2group <- data.frame(sample=character(0), group=character(0))
	for(g in groups$ID) {
	  members <- unlist(strsplit(groups[groups[,'ID']==g,'Members'],','))
	  sample2group <- rbind(sample2group, data.frame(sample=members, group=g))
	}
	
	# translate ensembl gene ids to HGNC symbols
	names(annotation)[names(annotation)=="Ensembl Gene ID"] <- "Ensembl"
	names(annotation)[names(annotation)=="Associated Gene Name"] <- "HGNC"
	
	hgncs <- get.parameter(cf, 'hgnc', 'string')
	
	genesViaParameter <- !is.null(hgncs) && !is.na(hgncs) && hgncs != ""
	if (genesViaParameter) {  # gene list (HGNC symbols) provided via parameter?
	  genes <- data.frame(HGNC=unlist(strsplit(hgncs, ',')))
	  genes <- merge(genes, annotation[,c("Ensembl", "HGNC")], all.x=T)
	} else {  # no: use genes (ensembl ids) provided via table input port
	  geneIds <- CSV.read(get.input(cf, "geneIds"))
	  genes <- merge(geneIds, annotation[,c("Ensembl", "HGNC")], all.x=T)
	  genes$HGNC[is.na(genes$HGNC) | genes$HGNC==""] <- genes$Ensembl
	}
	genes <- merge(genes, expr, by.x="Ensembl", by.y="row.names")
	
	# melt dataframe into shape suitable for plotting and add sample groups
	library(reshape)
	gexpr <- melt(genes, id.vars=c("Ensembl", "HGNC"))
	names(gexpr)[names(gexpr)=="variable"] <- "sample"
	gexpr <- merge(gexpr, sample2group, all.X=T)
	
	# before plotting, preserve input gene order in lattice output (note: factor levels determines panel order)
	# original order might have been destroyed by merge
	if (genesViaParameter) {
	  genes <- genes[match(unlist(strsplit(hgncs, ',')), genes$HGNC),] 
	} else {
	  genes <- genes[match(geneIds$Ensembl, genes$Ensembl),] # original order was destroyed by merge, so let's restore it first
	}                        	
	gexpr$HGNC <- factor(as.character(gexpr$HGNC), levels=genes$HGNC)
	
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
	
	# plot page by page
	library(lattice)
	
	nRow <- get.parameter(cf, 'nRow', 'int')
	nCol <- get.parameter(cf, 'nCol', 'int')
		
	splits <- split(levels(gexpr$HGNC), (0:(length(levels(gexpr$HGNC))-1)) %/% (nRow*nCol) + 1)

	width <- get.parameter(cf, 'width', 'float')
	height <- get.parameter(cf, 'height', 'float')
	
	rowheight <- height/nRow
	
	for(pageno in 1:length(splits)) {
	  gexpr.thispage <- gexpr[gexpr$HGNC %in% splits[[pageno]],]
	  nRow.thispage <- min(nRow, ((length(unique(gexpr.thispage$HGNC))-1) %/% nCol) + 1)
	  
	  # get output file name
	  plot.file <- sprintf('geneboxplot-%s-figure%d.pdf', instance.name, pageno)
	  
	  pdf(file.path(out.dir, plot.file), height=rowheight*nRow.thispage, width=width)
	  print(bwplot(value~group | HGNC, data=gexpr.thispage, 
	               layout=c(nCol,nRow.thispage),
	               par.strip.text=list(cex=0.5),
	               notch=FALSE,
	               as.table=TRUE,
	               ylab="DESeq2 normalized expression",
	               scales=list(x=list(rot=90)),
	               par.settings = list(box.umbrella=list(col="black"), box.rectangle = list(col="black")), 
	               panel=function(x,y,...){
	                 panel.grid()
	                 panel.bwplot(x,y,pch="|",do.out=FALSE, ...)
	                 panel.stripplot(x,y,jitter.data=TRUE,factor=0.8,pch=19,cex=0.2,...)
	               }))
	  dev.off()
	  
	  # generate latex string
	  caption <- get.parameter(cf, 'caption', 'string')
	  if (length(splits) > 1) {
	    caption <- paste0(caption, " Part ", pageno, " of ", length(splits), ".")
	  }
	  tex <- c(tex, latex.figure(plot.file, caption=caption))
	}
	
	latex.write.main(cf, 'document', tex)

	return(0)
}

main(execute)

