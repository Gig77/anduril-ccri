library(componentSkeleton)

execute <- function(cf) {

	# debug
	#rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/ikaros/results/anduril/execute/qcReport-gcbias/_command")

	# inputs
	
	countMatrix <- Matrix.read(get.input(cf, 'countMatrix'))
	sampleGroups <- CSV.read(get.input(cf, "sampleGroups"))
	annotation <- as.data.frame(Matrix.read(get.input(cf, "annotation")))
	
	# parameters 
	
	instance.name <- get.metadata(cf, 'instanceName')
	section.title <- get.parameter(cf, 'sectionTitle')
	section.type <- get.parameter(cf, 'sectionType')
	width <- get.parameter(cf, 'width', 'float')
	height <- get.parameter(cf, 'height', 'float')
	caption <- get.parameter(cf, 'caption', 'string')
	title <- get.parameter(cf, 'title', 'string')
	ymax <- get.parameter(cf, 'ymax', 'float')
	ymin <- get.parameter(cf, 'ymin', 'float')
	
	
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
	
	rownames(sample2group) <- sample2group$sample
	sample2group <- sample2group[sample2group$sample %in% colnames(countMatrix),]
	countMatrix <- countMatrix[,rownames(sample2group)]
	
	library(EDASeq)
	expressed <- apply(countMatrix, 1, function(x) sum(x) > 10)
	common <- intersect(rownames(countMatrix[expressed,]), rownames(annotation))
	data <- newSeqExpressionSet(
			counts=countMatrix[common,], 
			featureData = annotation[common,], 
			phenoData = data.frame(conditions=sample2group$group, row.names = sample2group$sample)
	)
	
	# prepare document
	tex <- character(0)
	tex <- c(tex, '\\clearpage')
	
	out.dir <- get.output(cf, 'document')
	dir.create(out.dir, recursive=TRUE)
	
	if (nchar(section.type) > 0) {
	  section.type=paste('\\',section.type,sep='') # if section type is defined, add escape in the beginning
	}
	if (nchar(section.title) > 0) {
	  tex <- c(tex, sprintf('%s{%s}\\label{%s}', section.type, section.title,	instance.name))
	}
	
	# plot
  plot.file <- sprintf('gcbias-%s.pdf', instance.name)
	pdf(file.path(out.dir, plot.file), height=height, width=width)
	library(RColorBrewer)
	palette(brewer.pal(8,"Paired"))
	biasPlot(data, "gc", log=TRUE, ylim=c(ymin,ymax), xlab="GC content (%)", ylab="Average read count", main=title)
	dev.off()
	
	# generate latex string
	tex <- c(tex, latex.figure(plot.file, caption=caption))
	latex.write.main(cf, 'document', tex)
	
	return(0)
}

main(execute)

