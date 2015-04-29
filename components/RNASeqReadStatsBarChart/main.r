library(componentSkeleton)
#library(gtools)

execute <- function(cf) {
	instance.name <- get.metadata(cf, 'instanceName')

	d <- read.delim(get.input(cf, 'stats'), row.names=1, check.names=F)
	d <- d[order(d$'uniquely mapped', decreasing=T),]
	
	document.dir <- get.output(cf, 'document')
	dir.create(document.dir, recursive=TRUE)
	plot.file <- sprintf('barchart-%s.pdf', instance.name)
	pdf(file.path(document.dir, plot.file))

	# all counts, grouped by sample; 8 sample per row; sorted by uniquely mapped reads in decreasing order
	samples.per.row <- get.parameter(cf, 'samplesPerRow', 'int')
	ymax <- get.parameter(cf, 'yMax', 'int')
	cols <- c("#332288", "#88CCEE", "#117733", "#DDCC77", "#CC6677", "#AA4499", "#44AA99", "#882255")
	
	chunks <- split(rownames(d), ceiling(seq_along(rownames(d)) / samples.per.row))
	par(mfrow=c(length(chunks), 1), mar=c(1,2,0.5,2), oma=c(1, 2, 2.5, 0))
	for(i in 1:length(chunks)) {
		chunk <- d[chunks[[i]],]
		
		# fill up missing samples in last row with empty plot
		for (j in nrow(chunk) : samples.per.row) {
			if (j == samples.per.row) { break }
			chunk <- rbind(chunk, c(0, 0, 0, 0, 0, 0, 0, 0))
			rownames(chunk)[j+1] <- paste0(rep(" ", j), collapse="")
		}
		
		barplot(t(chunk)/1000000, beside=T, xaxt='n', yaxt='n', ylim=c(0, ymax))
		abline(h=seq(0, ymax, by=5), col="gray", lty=2)
		par(cex.axis=0.6)
		bp <- barplot(t(chunk)/1000000, beside=T, col=cols, mgp=c(3,0,0), yaxt='n', ylim=c(0, ymax), add=T)
		box()
		#axis(side=1, at=bp, labels=rownames(chunk))
		axis(side=2, at=seq(0, ymax, by=1), labels=NA, tck=-0.01)
		axis(side=2, at=seq(0, ymax, by=5), labels=seq(0, ymax, by=5), mgp=c(3,0.7,0), las=2)		
	}
	mtext("number of reads x 1,000,000", outer=TRUE, 2)
	
	# create invisible overlay plot to put legend on top into outer margin
	par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
	plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
	legend("top", bg="white", colnames(d), fill=cols, ncol=4, cex=0.7, xpd=TRUE, bty="n")
	
	dev.off()
	
	# Prepare the document
	tex <- character()
	section.title <- get.parameter(cf, 'sectionTitle')
	section.type <- get.parameter(cf, 'sectionType')
	caption <- get.parameter(cf,'caption','string')
	myName <- get.metadata(cf, 'instanceName')
	if (nchar(section.type) > 0) {
		section.type=paste('\\',section.type,sep='') # if section type is defined, add escape in the beginning
	}
	if (nchar(section.title) > 0) {
		tex <- c(tex, sprintf('%s{%s}\\label{%s}', section.type, section.title,	myName))
	}
	
	tex <- c(tex, latex.figure(plot.file, caption=caption))		    			    
	latex.write.main(cf, 'document', tex)
	
	return(0)
}

main(execute)

