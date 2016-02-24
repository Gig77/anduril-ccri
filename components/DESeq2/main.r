library(componentSkeleton)

execute <- function(cf) {

  # debug
  #rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/fikret/results/anduril/execute/deseq_DTCfemale_DTCmale/_command")
  #stop("HERE!")
  
  instance.name <- get.metadata(cf, 'instanceName')	

  # inputs
  
  countMatrix  <- Matrix.read(get.input(cf, "countMatrix"))
  samples      <- CSV.read(get.input(cf, 'samples'))
  sampleGroups <- CSV.read(get.input(cf, 'sampleGroups'))
  
	# params
  
  caseGroup    <- get.parameter(cf, 'caseGroup',       type = 'string')
  controlGroup <- get.parameter(cf, 'controlGroup',    type = 'string')
	otherGroups  <- get.parameter(cf, "otherGroups",	 type = 'string')
	colSuffix    <- get.parameter(cf, 'colSuffix',       type = 'string')
	minReplicatesForReplace <- get.parameter(cf, "minReplicatesForReplace", type = 'int')
	design       <- get.parameter(cf, 'design',          type = 'string')
	reducedDesign       <- get.parameter(cf, 'reducedDesign',          type = 'string')
	coefficient   <- get.parameter(cf, 'coefficient',      type = 'string')
	cooksCutoff  <- get.parameter(cf, 'cooksCutoff',     type = 'float')
	trim  <- get.parameter(cf, 'trim',     type = 'float')
	
	stopifnot(controlGroup %in% sampleGroups$ID)
	stopifnot(caseGroup %in% sampleGroups$ID)
	
	if (is.na(design) || is.null(design) || design == "") {
		design = "~ group"
	}

	# prepare sample table for DESeq
  rownames(samples) <- samples$Alias
	samples$group <- NA

	controlSamples <- unlist(strsplit(sampleGroups[sampleGroups[,'ID']==controlGroup,'Members'],','))
	stopifnot(length(controlSamples) > 0)
	samples$group[samples$Alias %in% controlSamples] <- controlGroup

  caseSamples <- unlist(strsplit(sampleGroups[sampleGroups[,'ID']==caseGroup,'Members'],','))
  stopifnot(length(caseSamples) > 0)
  samples$group[samples$Alias %in% caseSamples] <- caseGroup

  # check
  samplesNotFound <- c(caseSamples, controlSamples)[!c(caseSamples, controlSamples) %in% colnames(countMatrix)]
  if (length(samplesNotFound) > 0) {
    stop("ERROR: The following samples were not present in count matrix: ", paste(samplesNotFound, sep = ","))
  }
  
  otherGroupsSplit <- NULL
  if (otherGroups != "") {
    otherGroupsSplit <- unlist(strsplit(otherGroups, ','))
    for (otherGroup in otherGroupsSplit) {
      otherGroupMembers <- unlist(strsplit(sampleGroups[sampleGroups[,'ID']==otherGroup,'Members'],','))
      stopifnot(length(otherGroupMembers) > 0)
      samples$group[samples$Alias %in% otherGroupMembers] <- otherGroup
    }
  }  
  
  samples <- samples[!is.na(samples$group),]
  samples$group <- factor(samples$group, levels=c(controlGroup, caseGroup, otherGroupsSplit))
  samples <- samples[order(as.integer(samples$group)),]
  countMatrix <- countMatrix[,rownames(samples)]

	print(sprintf("Groups in model: %s", paste(levels(samples$group), collapse=",")))
  print(sprintf("Contrasts: %s,%s", caseGroup, controlGroup))
  
  # DESeq2
  
  library(DESeq2)
  
  # read in count matrix
  dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = samples, design = as.formula(design))
  
  # estimate size factors
  dds <- estimateSizeFactors(dds)
  counts.norm <- as.data.frame(counts(dds, normalized=T))
  
  # estimate dispersion
  dds <- estimateDispersions(dds)
  
  # fit model
  if (reducedDesign == "") {
	  print("Performing negative binomial Wald test")
	  dds <- nbinomWaldTest(dds)
  } else {
	  print(paste0("Performing negative binomial LRT test with reduced design '", reducedDesign, "'"))
	  dds <- nbinomLRT(dds, reduced=as.formula(reducedDesign))
  }
  
  # replace outliers and fit again
  # NOTE: we have to do it this way because 'cooksCutoff' parameter cannot be specified with DESeq function
  # (it always defaults to 0.99)
  # original call instead of the above: dds <- DESeq(dds, minReplicatesForReplace=ifelse(minReplicatesForReplace > 0, minReplicatesForReplace, Inf))
  
  if (minReplicatesForReplace > 0) {
    print(paste0("Replacing outliers and re-fitting model: minReplicatesForReplace=", minReplicatesForReplace, ", cooksCutoff=", cooksCutoff))
    dds <- replaceOutliers(dds, trim = trim, cooksCutoff = cooksCutoff, minReplicates = minReplicatesForReplace)
	if (reducedDesign == "") {
		dds <- nbinomWaldTest(dds)
	} else {
		dds <- nbinomLRT(dds, reduced=as.formula(reducedDesign))
	}
}
  
  # extract results using case and control group as contrasts
	if (coefficient == "") {
		res <- results(dds, contrast=c("group", caseGroup, controlGroup), cooksCutoff=FALSE)
	} else {
		res <- results(dds, name=coefficient, cooksCutoff=FALSE)
	}
	results.out <- data.frame(ids=rownames(res), as.data.frame(res))
	
	# add normalized sample counts to output and compute separate means for experiment and control group	
	results.out <- merge(results.out, subset(counts.norm, select=colnames(countMatrix)), by.x="ids", by.y="row.names", all.x=T)
	results.out$baseMean <- rowMeans(results.out[,c(caseSamples, controlSamples)], na.rm=T)
	if (length(caseSamples) > 1) { results.out$baseMeanE <- rowMeans(results.out[,caseSamples], na.rm=T) } else { results.out$baseMeanE <-results.out[,caseSamples] }
	if (length(controlSamples) > 1) { results.out$baseMeanC <- rowMeans(results.out[,controlSamples], na.rm=T) } else { results.out$baseMeanC <-results.out[,controlSamples] }
	
	coln1 <- c('log2FoldChange', 'pvalue', 'padj', 'baseMean', 'baseMeanE', 'baseMeanC', 'lfcSE', 'stat')
	coln2 <- c('fc',             'p',      'q',    'meanExpr', 'meanExprE', 'meanExprC', 'fcSE',  'stat')
	coln2 <- paste(coln2, colSuffix, sep='')
	colnames(results.out)[match(coln1,colnames(results.out))] <- coln2
	results.out <- results.out[,c("ids", coln2, colnames(countMatrix))]
	results.out <- results.out[order(results.out[,paste0("q", colSuffix)]),]

	CSV.write(get.output(cf, 'results'), results.out)
	
	return(0)
}

main(execute)

