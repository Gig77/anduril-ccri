library(componentSkeleton)

execute <- function(cf) {
  
  # debug
  #rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/fikret/results/anduril/execute/keggDTCtum_vs_TUMorange-kegg/_command")
  #stop("HERE!")
  
  # metadata
  instance.name <- get.metadata(cf, 'instanceName')	
  
  # inputs
  
  degs <- CSV.read(get.input(cf, 'degs'))
  universe <- CSV.read(get.input(cf, 'universe'))
  
  # outputs
  report.dir <- get.output(cf, 'document')
  pathways.out <- get.output(cf, 'pathways')
  
  # params
  fdrCutoff <- get.parameter(cf, 'fdrCutoff', type = 'float')
  topN <- get.parameter(cf, 'topN', type = 'int')
  fcColumnName <- get.parameter(cf, 'fcColumnName', type = 'string')
  caption <- get.parameter(cf, 'caption',    type = 'string')
  section.title <- get.parameter(cf, 'sectionTitle', type = 'string')
  section.type <- get.parameter(cf, 'sectionType', type = 'string')
  
  # prepare Latex document
  dir.create(report.dir, recursive=TRUE)
  tex <- character()
  #tex <- c(tex, '\\clearpage')
  if (nchar(section.type) > 0) {
    section.type=paste('\\',section.type,sep='') # if section type is defined, add escape in the beginning
  }
  if (nchar(section.title) > 0) {
    tex <- c(tex, sprintf('%s{%s}\\label{%s}', section.type, section.title,	instance.name))
  }
  
  empty.df <- data.frame(ID=character(0), Description=character(0), GeneRatio=character(0), BgRatio=character(0), pvalue=numeric(0), p.adjust=numeric(0), qvalue=numeric(0), Count=integer(0), dir=character(0), geneID=character(0))
  
  # not genes left for enrichment? --> write empty output
  if (nrow(degs) <= 1) {
    latex.write.main(cf, 'document', tex)
    CSV.write(pathways.out, empty.df, first.cell = "ID")
    return(0)
  }

  # map Ensembl to Entrez ids
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  degs$entrez <- mapIds(org.Hs.eg.db, keys=as.character(degs[,1]), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  degs <- degs[!is.na(degs$entrez) & degs$entrez != "",]
  universe$entrez <- mapIds(org.Hs.eg.db, keys=as.character(universe[,1]), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  universe <- universe[!is.na(universe$entrez) & universe$entrez != "",]
  
  # enrichment
  library(clusterProfiler)
  
  # up-regulated genes
  entrez.up <- degs$entrez[degs[,fcColumnName] > 0]
  if (length(entrez.up) > 1) {
    ek <- enrichKEGG(entrez.up, organism="human", pvalueCutoff=1, pAdjustMethod="BH", universe=universe$entrez, qvalueCutoff = fdrCutoff)
    if (!is.null(ek) && !is.na(ek)) {
      enriched.up <- ek@result
      if (nrow(enriched.up) > 0) {
        enriched.up$dir <- "up"
      } else {
        enriched.up <- cbind(enriched.up, data.frame(dir=character(0)))
      }
    } else {
      enriched.up <- empty.df
    }
  } else {
    enriched.up <- empty.df
  }

  # down-regulated genes
  entrez.dn <- degs$entrez[degs[,fcColumnName] < 0]
  if (length(entrez.dn) > 1) {
    ek <- enrichKEGG(entrez.dn, organism="human", pvalueCutoff=1, pAdjustMethod="BH", universe=universe$entrez, qvalueCutoff = fdrCutoff)
    if (!is.null(ek) && !is.na(ek)) {
      enriched.dn <- ek@result
      if (nrow(enriched.dn) > 0) {
        enriched.dn$dir <- "down"
      } else {
        enriched.dn <- cbind(enriched.dn, data.frame(dir=character(0)))
      }
    } else {
      enriched.dn <- empty.df
    }
  } else {
    enriched.dn <- empty.df
  }
  
  # up- and down-regulated genes combined
  ek <- enrichKEGG(degs$entrez, organism="human", pvalueCutoff=1, pAdjustMethod="BH", universe=universe$entrez, qvalueCutoff = fdrCutoff)
  if (!is.null(ek) && !is.na(ek)) {
    enriched.both <- ek@result
    if (nrow(enriched.both) > 0) {
      enriched.both$dir <- "combined"
    } else {
      enriched.both <- cbind(enriched.both, data.frame(dir=character(0)))
    }
  } else {
    enriched.both <- empty.df
  }
  
  # combined directional and undirectional enrichments, sort, de-duplicate, re-order columns
  enriched <- rbind(enriched.up, enriched.dn, enriched.both)
  enriched <- enriched[order(enriched$pvalue),]
  enriched <- enriched[!duplicated(enriched$ID),] 
  enriched <- enriched[,c(colnames(enriched)[colnames(enriched)!="geneID"], "geneID")]
  
  # filter for signaling and metabolic pathways
  library("gage")
  data(kegg.gs)
  kegg.sets.hs <- kegg.gsets(species="hsa", id.type="kegg")
  kegg.sets.hs.sigmet <- gsub(" .*", "", names(kegg.sets.hs$kg.sets[kegg.sets.hs$sigmet.idx]))
  enriched <- enriched[enriched$ID %in% kegg.sets.hs.sigmet,]
  
  # map Entrez back to Ensembl ids
  if (nrow(enriched) > 0) {
    enriched$geneID <- sapply(strsplit(as.character(enriched$geneID), "/"), function(x) paste0(degs$ids[degs$entrez %in% x], collapse=","))
  }
  
  CSV.write(pathways.out, enriched, first.cell = "ID")
  
  # draw pathway diagrams
  library(pathview)
  foldchanges <- degs[,fcColumnName]
  names(foldchanges) <- degs$entrez
  foldchanges <- foldchanges[!is.na(degs$entrez)]
  if (nrow(enriched) > 0) {
    for (pw in enriched$ID[1:min(length(enriched$ID), topN)]) {
      print(paste("Plotting pathway", pw))
      pw.id <- gsub("hsa", "", pw)
      pw.file <- paste0(pw, ".pathview.png")
      plot.file <- sprintf('%s-%s.png', instance.name, pw)
      
      setwd(dirname(report.dir))
      pv.out <- pathview(gene.data = foldchanges, pathway.id = pw.id, species = "hsa", low=list(gene="#6666FF"), mid=list(gene="#CCCCCC"), high=list(gene="#FF0000"), kegg.dir = dirname(report.dir), same.layer=FALSE)
      system(paste0("mv ", getwd(), "/", pw.file, " ", report.dir, "/", plot.file))
      if (!file.exists(paste0(report.dir, "/", plot.file))) {
        print(paste0("WARNING: Could not create pathway diagram for '", pw, "'. Skipped!"))
        next
      }
      
      caption.prefix <- paste0("KEGG pathway '", pw, " ", enriched$Description[enriched$ID==pw], "'.")
      if (enriched$dir[enriched$ID==pw] == "up") {
        caption.suffix <- sprintf("Highest enrichment (q=%.2g) was found in up-regulated genes.", enriched$qvalue[enriched$ID==pw])
      } else if (enriched$dir[enriched$ID==pw] == "down") {
        caption.suffix <- sprintf("Highest enrichment (q=%.2g) was found in down-regulated genes.", enriched$qvalue[enriched$ID==pw])
      } else {
        caption.suffix <- sprintf("Highest enrichment (q=%.2g) was found considering both up- and down-regulated genes.", enriched$qvalue[enriched$ID==pw])
      }
      caption.suffix <- paste(caption.suffix, "Up- and down-regulated genes are colored in red and blue, respectively.")
      caption.combined <- ifelse(caption != "", paste(caption.prefix, caption, caption.suffix), paste(caption.prefix, caption.suffix))
      tex <- c(tex, latex.figure(plot.file, caption=caption.combined))
    }
  } else {
    tex <- c(tex, sprintf("No significantly enriched KEGG pathways with FDR cut-off %.2g found.", fdrCutoff))
  }
  
  latex.write.main(cf, 'document', tex)
  
  return(0)
}

main(execute)

