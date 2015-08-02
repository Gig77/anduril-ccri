library(componentSkeleton)

execute <- function(cf) {

  # debug
  #rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/iamp/results/anduril/execute/GSEAHeatmap/case1/component/_command")
  #rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/iamp/results/anduril/execute/gseaReportMSigDB-heatmapPositional/_command")

  instance.name <- get.metadata(cf, 'instanceName')	
	
	# Params
  caption <- get.parameter(cf, 'caption',    type = 'string')
  section.title <- get.parameter(cf, 'sectionTitle', type = 'string')
  section.type <- get.parameter(cf, 'sectionType', type = 'string')
  regexGeneSetNames <- get.parameter(cf, 'regexGeneSetNames',    type = 'string')
  regexGeneSetNamesExclude <- get.parameter(cf, 'regexGeneSetNamesExclude',    type = 'string')
  regexCategoryNames <- get.parameter(cf, 'regexCategoryNames',    type = 'string')
  regexCategoryNamesExclude <- get.parameter(cf, 'regexCategoryNamesExclude',    type = 'string')
  sigCutoff <- get.parameter(cf, 'sigCutoff',    type = 'float')
  hsigCutoff <- get.parameter(cf, 'hsigCutoff',    type = 'float')
  nesCutoff <- get.parameter(cf, 'NESCutoff', type = 'float')
  filterExpr <- get.parameter(cf, 'RFilterExpr', type = 'string')
  cexCol <- get.parameter(cf, 'cexCol', type = 'float')
  cexRow <- get.parameter(cf, 'cexRow', type = 'float')
  title <- get.parameter(cf, 'title',    type = 'string')
  figureHeight <- get.parameter(cf, 'figureHeight',    type = 'float')
  pdfHeight <- get.parameter(cf, 'pdfHeight',    type = 'float')
  pdfWidth <- get.parameter(cf, 'pdfWidth',    type = 'float')
  
	# Inputs
	enrichedUp  <- Array.read(cf,"enrichedUp")
	enrichedDown  <- Array.read(cf,"enrichedDown")
	
	stopifnot(nrow(enrichedUp)==nrow(enrichedDown))
	stopifnot(nrow(enrichedUp) > 0)
	stopifnot(nrow(enrichedDown) > 0)
	
	# read and merge gene sets
	merged <- read.delim(enrichedUp[1,"File"], colClasses = c("character", "character", "character", "NULL", "NULL", "NULL", "numeric", "NULL", "numeric", "NULL", "NULL", "NULL", "NULL"))
	names(merged)[4:ncol(merged)] <- paste0(names(merged)[4:ncol(merged)], ".", enrichedUp[1,"Key"], ".up")
	if (nrow(enrichedUp) > 1) {
	  for (i in 2:nrow(enrichedUp)) {
	    d <- read.delim(enrichedUp[i,"File"], colClasses = c("character", "character", "character", "NULL", "NULL", "NULL", "numeric", "NULL", "numeric", "NULL", "NULL", "NULL", "NULL"))
	    names(d)[4:ncol(d)] <- paste0(names(d)[4:ncol(d)], ".", enrichedUp[i,"Key"], ".up")
	    merged <- merge(merged, d, by=1:3, all=T)
	  }
	}
	for (i in 1:nrow(enrichedDown)) {
	  d <- read.delim(enrichedDown[i,"File"], colClasses = c("character", "character", "character", "NULL", "NULL", "NULL", "numeric", "NULL", "numeric", "NULL", "NULL", "NULL", "NULL"))
	  names(d)[4:ncol(d)] <- paste0(names(d)[4:ncol(d)], ".", enrichedDown[i,"Key"], ".down")
	  merged <- merge(merged, d, by=1:3, all=T)
	}
	rownames(merged) <- merged$NAME
	
	# for each up/down pair, determine higher (squared) NES and smaller q-val 
	for (cmp in enrichedUp$Key) {
	  merged[,cmp] <- ifelse(!is.na(merged[,paste0("NES.", cmp, ".up")]), merged[,paste0("NES.", cmp, ".up")], merged[,paste0("NES.", cmp, ".down")])
	  merged[,cmp][is.na(merged[,cmp])] <- 0
	  merged[,cmp] <- ifelse(merged[,cmp] > 0, merged[,cmp]*merged[,cmp], -merged[,cmp]*merged[,cmp]) # square NES to increase color contrast
	  merged[,paste0("FDR.q.val.", cmp)] <- pmin(merged[,paste0("FDR.q.val.", cmp, ".up")], merged[,paste0("FDR.q.val.", cmp, ".down")], na.rm=T)
	  merged[,paste0("FDR.q.val.", cmp)][is.na(merged[,paste0("FDR.q.val.", cmp)])] <- 1
	}

	# determine best FDR and NES across all comparisons
	merged$FDR.best <- apply(merged[,grepl("FDR.q.val", names(merged))], 1, min, na.rm=TRUE)
	merged$NES.best <- apply(merged[,grepl("^NES.", names(merged))], 1, function (x) { max(abs(x), na.rm=TRUE)} )
	
	# geneset filtering
	genesets <- merged
	if (!is.na(regexCategoryNames) && !is.null(regexCategoryNames) && nchar(regexCategoryNames) > 0) {
	  genesets <- genesets[grepl(regexCategoryNames, genesets$CATEGORY, perl=T),]
	}
	if (!is.na(regexCategoryNamesExclude) && !is.null(regexCategoryNamesExclude) && nchar(regexCategoryNamesExclude) > 0) {
	  genesets <- genesets[!grepl(regexCategoryNamesExclude, genesets$CATEGORY, perl=T),]
	}
	if (!is.na(regexGeneSetNames) && !is.null(regexGeneSetNames) && nchar(regexGeneSetNames) > 0) {
	  genesets <- genesets[grepl(regexGeneSetNames, genesets$NAME, perl=T),]
	}
	if (!is.na(regexGeneSetNamesExclude) && !is.null(regexGeneSetNamesExclude) && nchar(regexGeneSetNamesExclude) > 0) {
	  genesets <- genesets[!grepl(regexGeneSetNamesExclude, genesets$NAME, perl=T),]
	}
	
	if (nchar(filterExpr) > 0) {
	  genesets <- with(genesets, genesets[eval(parse(text=filterExpr)),])
	} else {
    genesets <- genesets[genesets$FDR.best <= sigCutoff,]
  }
	
	if (nesCutoff > 0) {
	  genesets <- genesets[genesets$NES.best >= nesCutoff,]
	}

	# prepare output directory
	report.dir <- get.output(cf, 'document')
	dir.create(report.dir, recursive=TRUE)
	
	# prepare Latex document
	tex <- character()
	tex <- c(tex, '\\clearpage')
	if (nchar(section.type) > 0) {
	  section.type=paste('\\',section.type,sep='') # if section type is defined, add escape in the beginning
	}
	if (nchar(section.title) > 0) {
	  tex <- c(tex, sprintf('%s{%s}\\label{%s}', section.type, section.title,	instance.name))
	}

	if (nrow(genesets) < 2) {
	  print("WARNING: Not enough (< 2) gene sets to produce cluster heatmap. Skipped.")
	  tex <- c(tex, "Not enough enriched gene sets in this category. No heatmap produced.")
	} else {
	  # plot heatmap	
	  library("RColorBrewer")
	  library("gplots")
	  
	  hmdata <- genesets[,enrichedUp$Key]
	  
	  hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
	  fdrs <- as.matrix(merged[rownames(hmdata), grep("(up|down)$", grep("^FDR.q.val", names(merged), value=T), value=T, invert=T)])
	  sig <- fdrs ; sig[,] <- NA
	  sig[fdrs <= sigCutoff & fdrs > hsigCutoff] <- "*"
	  sig[fdrs <= hsigCutoff] <- "**"
	  
	  plot.file <- sprintf('%s-gseaheatmap.pdf', instance.name)
	  pdf(file.path(report.dir, plot.file), height=pdfHeight, width=pdfWidth)
	  heatmap.2(as.matrix(hmdata), Colv=F, Rowv=T, dendrogram="row", trace="none", col=rev(hmcol), margin=c(10, 25), cexCol=cexCol, cexRow=cexRow, keysize=0.7, 
	            colsep=seq(1:ncol(hmdata)), rowsep=seq(1:nrow(hmdata)), sepcolor="grey92", sepwidth=c(0.005,0.005),
	            cellnote=sig, notecol='white',
	            main=sprintf("%s\n(* FDR <= %.4g, ** FDR <= %.4g%s)", title, sigCutoff, hsigCutoff, ifelse(nesCutoff > 0, sprintf(",\nabs(NES) >= %.1f", nesCutoff), "")),
	            key.title="NES^2", key.xlab="", key.ylab="")
	  dev.off()
	  
	  tex <- c(tex, latex.figure(plot.file, caption=caption, image.height=figureHeight))  
	}
	
	latex.write.main(cf, 'document', tex)
	
	return(0)
}

main(execute)

