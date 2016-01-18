library(componentSkeleton)

execute <- function(cf) {

  # debug
  #rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/iamp/results/anduril/execute/rnk_iAMPvsNoniAMP/_command")
  #stop("HERE!")
  
  instance.name <- get.metadata(cf, 'instanceName')	

	deg <- CSV.read(get.input(cf, 'deg'))
	annotation <- CSV.read(get.input(cf, 'annotation'))

	col.p <- get.parameter(cf, "colP")
	col.fc <- get.parameter(cf, "colFC")
	min.p <- get.parameter(cf, "minP", "float")
	rank.by <- get.parameter(cf, "rankBy", "string")
	ignore.direction <- get.parameter(cf, "ignoreDirection", "boolean")
	
	# map IDs to name
	deg.named <- merge(deg[,c("ids", col.p, col.fc)], annotation[,c(1,2)], by.x="ids", by.y=names(annotation)[1])

	# filtering
	deg.named <- deg.named[!is.na(deg.named[,col.p]) & !is.na(deg.named[,col.fc]) & deg.named[,4] != "",]
	deg.named <- deg.named[deg.named[,col.p] <= min.p,]
	
	# de-duplication: if multiple gene IDs mapped to same name, keep only entry with smallest p
	deg.named <- deg.named[ave(deg.named[,col.p], deg.named[,4], FUN=min) == deg.named[,col.p],]

	# convert p-values to scores where sign indicates directionality (up or downregulated)
	if (rank.by == "p") {
		if (ignore.direction) {
			deg.named$score <- -log(deg.named[,col.p],10)
		} else {
			deg.named$score <- ifelse(deg.named[,col.fc]>0, -log(deg.named[,col.p],10), log(deg.named[,col.p],10))
		}
	} else if (rank.by == "fc") {
		if (ignore.direction) {
			deg.named$score <- abs(deg.named[,col.fc])
		} else {
			deg.named$score <- deg.named[,col.fc]
		}
	} else if (rank.by == "pi") {
		if (ignore.direction) {
			deg.named$score <- abs(deg.named[,col.fc]) * -log(deg.named[,col.p],10)
		} else {
			deg.named$score <- deg.named[,col.fc] * -log(deg.named[,col.p],10)
		}
	} else {
		stop(sprintf("ERROR: Unknown value for parameter 'rankBy': %s", rank.by))
	}
	
	# keep only gene name and score; uppercase gene name to allow GSEA with mouse HGNC symbols
	rnk <- deg.named[order(deg.named$score,decreasing=T), c(4, 5)]
	rnk[,1] <- toupper(rnk[,1])
	
	# write output
	out.file <- get.output(cf, 'rnk')
	write.table(rnk, out.file, col.names=F, row.names=F, quote=F, sep="\t")

	return(0)
}

main(execute)

