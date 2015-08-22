library(componentSkeleton)

execute <- function(cf) {

  # debug
  #rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/iamp/results/anduril/execute/GSEAOverlap/case1/component/_command")
  #rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/helena_veronika/results/anduril/execute/gseaReportOverlap-dendrogramoeRHDvsEmptyDn/_command")
  
  instance.name <- get.metadata(cf, 'instanceName')	
	
  # Params
  caption <- get.parameter(cf, 'caption',    type = 'string')
  section.title <- get.parameter(cf, 'sectionTitle', type = 'string')
  section.type <- get.parameter(cf, 'sectionType', type = 'string')
  regexGeneSetNames <- get.parameter(cf, 'regexGeneSetNames',    type = 'string')
  regexGeneSetNamesExclude <- get.parameter(cf, 'regexGeneSetNamesExclude',    type = 'string')
  regexCategoryNames <- get.parameter(cf, 'regexCategoryNames',    type = 'string')
  regexCategoryNamesExclude <- get.parameter(cf, 'regexCategoryNamesExclude',    type = 'string')
  similarityCutoffTree <- get.parameter(cf, 'similarityCutoffTree',    type = 'float')
  similarityCutoffTable <- get.parameter(cf, 'similarityCutoffTable',    type = 'float')
  sigCutoff <- get.parameter(cf, 'sigCutoff',    type = 'float')
  hsigCutoff <- get.parameter(cf, 'hsigCutoff',    type = 'float')
  nesCutoff <- get.parameter(cf, 'nesCutoff',    type = 'float')
  topN <- get.parameter(cf, 'topN',    type = 'int')
  cexLabel <- get.parameter(cf, 'cexLabel', type = 'float')
  
  # Inputs
  enrichedGeneSets  <- Array.read(cf, "enrichedGeneSets")

	stopifnot(nrow(enrichedGeneSets) > 0)

	# read and merge gene sets
	sets <- data.frame(NAME=character(0),	
                     DESCRIPTION=character(0),	
                     CATEGORY=character(0),	
	                   LINKOUT=character(0),	
	                   SIZE=integer(0),	
	                   ES=numeric(0), 
	                   NES=numeric(0), 
	                   NOM.p.val=numeric(0), 
	                   FDR.q.val=numeric(0), 
	                   FWER.p.val=numeric(0), 
	                   RANK.AT.MAX=integer(0),
	                   LEADING.EDGE=character(0),
                     CORE_GENES=character(0),
	                   stringsAsFactors = FALSE)
  
	for (file in enrichedGeneSets$File) {
    d <- read.delim(file, stringsAsFactors=FALSE)
    sets <- rbind(sets, d)
  }

	# filter gene sets
  sets <- sets[sets$FDR.q.val <= sigCutoff,]
  if (nesCutoff > 0) {
	  sets <- sets[sets$NES >= nesCutoff,]
  }
  if (!is.na(regexCategoryNames) && !is.null(regexCategoryNames) && nchar(regexCategoryNames) > 0) {
    sets <- sets[grepl(regexCategoryNames, sets$CATEGORY, perl=T),]
  }
  if (!is.na(regexCategoryNamesExclude) && !is.null(regexCategoryNamesExclude) && nchar(regexCategoryNamesExclude) > 0) {
    sets <- sets[!grepl(regexCategoryNamesExclude, sets$CATEGORY, perl=T),]
  }
  if (!is.na(regexGeneSetNames) && !is.null(regexGeneSetNames) && nchar(regexGeneSetNames) > 0) {
    sets <- sets[grepl(regexGeneSetNames, sets$NAME, perl=T),]
  }
  if (!is.na(regexGeneSetNamesExclude) && !is.null(regexGeneSetNamesExclude) && nchar(regexGeneSetNamesExclude) > 0) {
    sets <- sets[!grepl(regexGeneSetNamesExclude, sets$NAME, perl=T),]
  }

  # prepare Latex document
  tex <- character()
  tex <- c(tex, '\\clearpage')
  if (nchar(section.type) > 0) {
    section.type=paste('\\',section.type,sep='') # if section type is defined, add escape in the beginning
  }
  if (nchar(section.title) > 0) {
    tex <- c(tex, sprintf('%s{%s}\\label{%s}', section.type, section.title,	instance.name))
  }

  overlaps <- data.frame(set1=character(0), set2=character(0), desc1=character(0), desc2=character(0), linkout1=character(0), linkout2=character(0), size1=integer(0), size2=integer(0), shared=integer(0), jacc=numeric(0), NES1=numeric(0), NES2=numeric(0), FDR1=numeric(0), FDR2=numeric(0), shared_core_genes=character(0), stringsAsFactors = F)
  
  if (nrow(sets) < 2) {
    print("WARNING: Not enough (<2) overlapping gene sets to produce dendrogram. Skipped.")
    tex <- c(tex, "No overlapping gene sets found.")
  } else {
    # set up distance matrix based on Jaccard indices
  	dist <- matrix(0, nrow=nrow(sets), ncol=nrow(sets), dimnames=list(sets$NAME, sets$NAME))
  	
  	# remember which genes were found in which gene sets (and its enrichment)
  	coreGenes <- data.frame(gene=character(), set=character(), fdr=numeric(), stringsAsFactors = 0)
  	
  	for (i in c(1:nrow(sets))) {
  		g1 <- unlist(strsplit(sets$CORE_GENES[i], ", "))
  		
  		for (j in c(1:nrow(sets))) {
  			g2 <- unlist(strsplit(sets$CORE_GENES[j], ", "))
  			shared <- intersect(g1, g2)
  			jacc = length(shared) / length(union(g1, g2))
  			dist[sets$NAME[i], sets$NAME[j]] <- 1-jacc
  			if (i < j && jacc >= similarityCutoffTable) {
  				overlaps <- rbind(overlaps, data.frame(set1=sets$NAME[i], set2=sets$NAME[j], desc1=sets$DESCRIPTION[i], desc2=sets$DESCRIPTION[j], linkout1=sets$LINKOUT[i], linkout2=sets$LINKOUT[j], size1=length(g1), size2=length(g2), shared=length(shared), jacc=jacc, NES1=sets$NES[i], NES2=sets$NES[j], FDR1=sets$FDR.q.val[i], FDR2=sets$FDR.q.val[j], shared_core_genes=paste(intersect(g1, g2), collapse=", "), stringsAsFactors = F))
  			}
  		}
  		
  		for (g in g1) {
  		  coreGenes[nrow(coreGenes)+1,"gene"] <- g ; coreGenes[nrow(coreGenes),"set"] <- sets$NAME[i] ; coreGenes[nrow(coreGenes),"fdr"] <- sets$FDR.q.val[i]
  		}
  	}
  	overlaps <- overlaps[order(overlaps$jacc, decreasing=T),]
  	
    # prune tree: exclude gene sets with insufficient overlap
    diag(dist) <- 1
    include <- apply(dist, 1, min) <= 1-similarityCutoffTree
    dist.pruned <- dist[include, include]
	
	# prunte to top N sets if requested
  caption.suffix = ""
	if (topN > 0 & nrow(dist.pruned) > topN) {
	  caption.suffix = paste(" Of ", nrow(dist.pruned), " overlapping gene sets meeting the filtering criteria, only the ", topN, " most enriched gene sets are shown.")
	  sets.pruned <- sets[sets$NAME %in% colnames(dist.pruned),]
		sets.pruned <- sets.pruned[order(abs(sets.pruned$NES), decreasing=T),]
		include <- sets.pruned$NAME[1:topN]
		dist.pruned <- dist.pruned[include, include]
		overlaps <- overlaps[overlaps$set1 %in% include & overlaps$set2 %in% include,]
	}
	
    # prepare output directory
    report.dir <- get.output(cf, 'dendrogram')
    dir.create(report.dir, recursive=TRUE)
    
  	if (nrow(dist.pruned) < 2) {
  	  print("WARNING: Not enough (<2) overlapping gene sets to produce dendrogram. Skipped.")
  	  tex <- c(tex, "Not enough overlapping gene sets. No dendrogram produced.")
  	} else {
  	  
  	  # set label expansion factor
  	  if (cexLabel==0) {
  	    cexLabel = 1 - nrow(dist.pruned) / 130;
  	    if (cexLabel < 0.2) {
  	      cexLabel = 0.2
  	    }
  	  }
  	  
  	  # plot dendrogram
  	  plot.file <- sprintf('%s-tree.pdf', instance.name)
  	  pdf(file.path(report.dir, plot.file))
  	  rd <- as.dendrogram(hclust(as.dist(dist.pruned)), hang=0.02)

  	  # color dendrogram leaves if FDR <= hsigCutoff
  	  if (hsigCutoff >= 0) {
  	    labelCol <- function(x) {
  	      if (is.leaf(x) && sets$FDR.q.val[sets$NAME==attr(x, "label")] <= hsigCutoff) {
  	        attr(x, "nodePar") <- list(lab.col="red")
  	      }
  	      return(x)
  	    }
  	    rd <- dendrapply(rd, labelCol)
  	  }
  	  
  	  par(cex=cexLabel, mar=c(4,1,1,10)) 
  	  plot(rd, horiz=T, axes = F, xlab="Shared core enrichment genes (Jaccard index)", xlim=c(1,0), cex.lab=2.5-1.5*cexLabel)
  	  #plot(hclust(as.dist(dist)))
  	  axis(1, at=seq(0,1,0.1), labels=1-seq(0,1,0.1))
  	  dev.off()
  
  	  tex <- c(tex, latex.figure(plot.file, caption=paste0(caption, caption.suffix)))  
  	}
  }
  
	latex.write.main(cf, 'dendrogram', tex)
	CSV.write(get.output(cf, 'overlaps'), overlaps)
	
	# write enriched core genes with their gene sets
	if (nrow(sets) < 2) {
	  coreGenes.agg <- data.frame(gene=NA, rank=NA, freq=NA, 'gene_sets(FDR)'=NA, check.names=F)
	} else {
	  library(plyr)
    coreGenes$set_fdr <- paste0(coreGenes$set, "(", sprintf("%.1g", coreGenes$fdr), ")")
    coreGenes.agg <- ddply(coreGenes[order(coreGenes$gene, coreGenes$fdr, coreGenes$set),], ~gene, summarise, freq = length(gene), 'gene_sets(FDR)' = paste(set_fdr, collapse = ", "))
    coreGenes.agg <- coreGenes.agg[order(coreGenes.agg$freq, decreasing=T),]
    names(coreGenes.agg)[names(coreGenes.agg) == "gene"] <- "gene_rank"
    coreGenes.agg$gene <- sapply(strsplit(coreGenes.agg$gene_rank, "\\(|\\)"), "[[", 1)
    coreGenes.agg$rank <- sapply(strsplit(coreGenes.agg$gene_rank, "\\(|\\)"), "[[", 2)
	}

	CSV.write(get.output(cf, 'coreGeneFrequency'), coreGenes.agg[,c("gene", "rank", "freq", "gene_sets(FDR)")])	
	
	return(0)
}

main(execute)

