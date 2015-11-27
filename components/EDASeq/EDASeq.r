library(componentSkeleton)

execute <- function(cf) {

	# debug
	#rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/ikaros/results/anduril/execute/edaseq/_command")

  instance.name <- get.metadata(cf, 'instanceName')
  
	# inputs
	
	countMatrix <- Matrix.read(get.input(cf, 'countMatrix'))
	annotation <- as.data.frame(Matrix.read(get.input(cf, "annotation")))
	
	# parameters 
	
	minReadCount <- get.parameter(cf, 'minReadCount', 'int')
	method <- get.parameter(cf, 'method', 'string')
	
	library(EDASeq)
	
	# construct EDASeq expression count matrix using only expressed genes
	expressed <- apply(countMatrix, 1, function(x) sum(x) >= minReadCount)
	common <- intersect(rownames(countMatrix[expressed,]), rownames(annotation))
	data <- newSeqExpressionSet(counts=countMatrix[common,], featureData = annotation[common,])
	
	# within-lane normalization
	dataNorm.withinLane <- withinLaneNormalization(data, "gc", which=method, offset=TRUE, round=FALSE)
	CSV.write(get.output(cf, 'normalizedCountsWithinLane'), round(normCounts(dataNorm.withinLane)), first.cell = "ids")
	CSV.write(get.output(cf, 'offsetsWithinLane'), offst(dataNorm.withinLane), first.cell = "ids")

	# between-lane normalization
	dataNorm.betweenLane <- betweenLaneNormalization(dataNorm.withinLane, which=method, offset=TRUE, round=FALSE)
	CSV.write(get.output(cf, 'normalizedCountsBetweenLane'), round(normCounts(dataNorm.betweenLane)), first.cell = "ids")
	CSV.write(get.output(cf, 'offsetsBetweenLane'), offst(dataNorm.betweenLane), first.cell = "ids")

	# DESeq normalization factors computed from offsets
	normFactors <- exp(-1 * offst(dataNorm.betweenLane))
	normFactors <- normFactors / exp(rowMeans(log(normFactors)))
	CSV.write(get.output(cf, 'DESeqNormFactors'), normFactors, first.cell = "ids")

	#pdf("/mnt/projects/ikaros/results/gc-correction.edaseq.pdf", width=10, height=10)
	#par(mfrow=c(5,5), mar=c(0,0,0,0))
	#for (i in 1:ncol(countMatrix)) {
	#  plot(countMatrix[common,i]+0.1, round(normCounts(dataNorm.betweenLane))[common,i]+0.1, log="xy", cex=0.1, xlim=c(10,10000), ylim=c(10,10000), xlab="", ylab="", xaxt='n', yaxt='n', col="red") 
	#  abline(0, 1, col="black", lwd=3)
	#  text(10,8000,paste0(colnames(countMatrix)[i]), adj=0)
	#}
	#dev.off()
	
	return(0)
}

main(execute)

