library(componentSkeleton)

execute <- function(cf) {

  # debug
  #rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/helena_veronika/results/anduril/execute/deseq_oeERvsEmpty/_command")
  #stop("HERE!")
  
  instance.name <- get.metadata(cf, 'instanceName')	

  # inputs
  
  countFiles   <- Array.read(cf,"counts")
  samples      <- CSV.read(get.input(cf, 'samples'))
  sampleGroups <- CSV.read(get.input(cf, 'sampleGroups'))
  
	# params
  
	controlGroup <- get.parameter(cf, 'controlGroup',    type = 'string')
	caseGroup    <- get.parameter(cf, 'caseGroup',       type = 'string')
	otherGroups  <- get.parameter(cf, "otherGroups", type = 'string')
	label        <- get.parameter(cf, 'label',          type = 'string')
	minReplicatesForReplace <- get.parameter(cf, "minReplicatesForReplace", type = 'int')
	design       <- get.parameter(cf, 'design',          type = 'string')
	
	stopifnot(controlGroup %in% sampleGroups$ID)
	stopifnot(caseGroup %in% sampleGroups$ID)
	
	if (is.na(design) || is.null(design) || design == "") {
		design = "~ group"
	}

	# prepare sample table for DESeq

	sampleTable <- merge(samples, countFiles, by.x="Alias", by.y="Key", all.x = T)
	names(sampleTable)[ncol(sampleTable)] <- "fileName"
	sampleTable$group <- NA
  
	controlSamples <- unlist(strsplit(sampleGroups[sampleGroups[,'ID']==controlGroup,'Members'],','))
	stopifnot(length(controlSamples) > 0)
	sampleTable$group[sampleTable$Alias %in% controlSamples] <- controlGroup

  caseSamples <- unlist(strsplit(sampleGroups[sampleGroups[,'ID']==caseGroup,'Members'],','))
  stopifnot(length(caseSamples) > 0)
  sampleTable$group[sampleTable$Alias %in% caseSamples] <- caseGroup

  otherGroupsSplit <- NULL
  if (otherGroups != "") {
    otherGroupsSplit <- unlist(strsplit(otherGroups, ','))
    for (otherGroup in otherGroupsSplit) {
      otherGroupMembers <- unlist(strsplit(sampleGroups[sampleGroups[,'ID']==otherGroup,'Members'],','))
      stopifnot(length(otherGroupMembers) > 0)
      sampleTable$group[sampleTable$Alias %in% otherGroupMembers] <- otherGroup
    }
  }  
  
  sampleTable <- sampleTable[!is.na(sampleTable$group),]
  sampleTable$group <- factor(sampleTable$group, levels=c(controlGroup, caseGroup, otherGroupsSplit))
  sampleTable <- sampleTable[,c(c("Alias", "fileName", "group"), names(sampleTable)[!names(sampleTable) %in% c("Alias", "fileName", "group")])]


	print(sprintf("Groups in model: %s", paste(levels(sampleTable$group), collapse=",")))
  print(sprintf("Contrasts: %s,%s", caseGroup, controlGroup))
  
  # DESeq2
  
	library(DESeq2)
	ddsHTSeq  <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = "/", design = as.formula(design))
	ddsHTSeq  <- DESeq(ddsHTSeq, minReplicatesForReplace=ifelse(minReplicatesForReplace > 0, minReplicatesForReplace, Inf))
	res       <- results(ddsHTSeq, contrast=c("group", caseGroup, controlGroup), cooksCutoff=FALSE)
	results.out <- data.frame(ids=rownames(res), as.data.frame(res))
	
	# add normalized sample counts to output and compute separate means for experiment and control group	
	counts.norm <- as.data.frame(counts(ddsHTSeq, normalized=T))
	results.out <- merge(results.out, subset(counts.norm, select=c(caseSamples, controlSamples)), by.x="ids", by.y="row.names", all.x=T)
	results.out$baseMean <- rowMeans(results.out[,c(caseSamples, controlSamples)], na.rm=T)
	if (length(caseSamples) > 1) { results.out$baseMeanE <- rowMeans(results.out[,caseSamples], na.rm=T) } else { results.out$baseMeanE <-results.out[,caseSamples] }
	if (length(controlSamples) > 1) { results.out$baseMeanC <- rowMeans(results.out[,controlSamples], na.rm=T) } else { results.out$baseMeanC <-results.out[,controlSamples] }
	
	coln1 <- c('log2FoldChange', 'pvalue', 'padj', 'baseMean', 'baseMeanE', 'baseMeanC', 'lfcSE', 'stat')
	coln2 <- c('fc',             'p',      'q',    'meanExpr', 'meanExprE', 'meanExprC', 'fcSE',  'stat')
	coln2 <- paste(coln2, label, sep='')
	colnames(results.out)[match(coln1,colnames(results.out))] <- coln2
	results.out <- results.out[,c("ids", coln2, caseSamples, controlSamples)]
	results.out <- results.out[order(results.out[,paste0("q", label)]),]

	CSV.write(get.output(cf, 'results'), results.out)
	
	return(0)
}

main(execute)

