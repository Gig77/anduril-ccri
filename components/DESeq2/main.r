library(componentSkeleton)

execute <- function(cf) {

	# Meta data
	instance.name <- get.metadata(cf, 'instanceName')	
	
	# Params
	nameControl <- get.parameter(cf, 'nameControl',    type = 'string')
	nameCase    <- get.parameter(cf, 'nameCase',       type = 'string')
	label       <- get.parameter(cf, 'label',          type = 'string')

	# Inputs
	countFiles  <- Array.read(cf,"countFiles")
	samples     <- CSV.read(get.input(cf, 'samples'))
	
	sNames      <- countFiles$Key
	g1N         <- unlist(strsplit(samples[samples[,'ID']==nameControl,'Members'],','))
	g2N         <- unlist(strsplit(samples[samples[,'ID']==nameCase,'Members'],','))
	condition   <- rep(NA, length(sNames))
	condition[match(g1N,sNames)] <- nameControl
	condition[match(g2N,sNames)] <- nameCase
	condition   <- factor(x=condition, levels=c(nameControl,nameCase), ordered=FALSE)
	sampleTable <- data.frame(sampleName = sNames, fileName = countFiles$File, condition  = condition)
	sampleTable <- sampleTable[!is.na(condition),]	
	
	library(DESeq2)
	ddsHTSeq  <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "/", design = ~ condition)
	ddsHTSeq  <- DESeq(ddsHTSeq)
	res       <- results(ddsHTSeq, cooksCutoff=FALSE)
	results.out <- data.frame(ids=rownames(res), as.data.frame(res))
	
	# add normalized sample counts to output and compute separate means for experiment and control group	
	counts.norm <- as.data.frame(counts(ddsHTSeq, normalized=T))
	results.out <- merge(results.out, counts.norm, by.x="ids", by.y="row.names", all.x=T)
	if (length(g2N) > 1) { results.out$baseMeanE <- rowMeans(results.out[,g2N], na.rm=T) } else { results.out$baseMeanE <-results.out[,g2N] }
	if (length(g1N) > 1) { results.out$baseMeanC <- rowMeans(results.out[,g1N], na.rm=T) } else { results.out$baseMeanC <-results.out[,g1N] }
	
	coln1 <- c('pvalue', 'padj', 'baseMean', 'baseMeanE', 'baseMeanC', 'log2FoldChange', 'lfcSE')
	coln2 <- c('p',      'q',    'meanExpr', 'meanExprE', 'meanExprC', 'fc',             'fcSE')
	coln2 <- paste(coln2, label, sep='')
	colnames(results.out)[match(coln1,colnames(results.out))] <- coln2

	CSV.write(get.output(cf, 'results'), results.out)
	
	return(0)
}

main(execute)

