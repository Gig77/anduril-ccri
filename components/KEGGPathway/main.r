library(componentSkeleton)

execute <- function(cf) {

  # debug
  #rm(list=ls()) ; cf <- parse.command.file("/mnt/projects/ikaros/results/anduril/execute/keggIKM_vs_IKC/_command")
  #stop("HERE!")

  # metadata
  instance.name <- get.metadata(cf, 'instanceName')	

  # inputs
  
  degs <- CSV.read(get.input(cf, 'degs'))

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
  
  library("AnnotationDbi")
  library("org.Hs.eg.db")
  
  # get fold-changes
  degs$entrez <- mapIds(org.Hs.eg.db, keys=as.character(degs[,1]), column="ENTREZID", keytype="ENSEMBL", multiVals="first") 
  foldchanges <- degs[,fcColumnName]
  names(foldchanges) <- degs$entrez
  foldchanges <- foldchanges[!is.na(names(foldchanges))]

  # enrichment
  library("gage")
  data(kegg.gs)
  kegg.sets.hs <- kegg.gsets(species="hsa", id.type="kegg")
  keggres <- gage(foldchanges, gsets=kegg.sets.hs$kg.sets[kegg.sets.hs$sigmet.idx], same.dir=FALSE, ref = NULL, samp = NULL)
  keggres.samedir <- gage(foldchanges, gsets=kegg.sets.hs$kg.sets[kegg.sets.hs$sigmet.idx], same.dir=TRUE, ref = NULL, samp = NULL)
  
  # combine three results, take most significant one
  keggres.df <- as.data.frame(keggres$greater) ; keggres.df$dir <- "combined" ; keggres.df$name <- rownames(keggres.df) ; rownames(keggres.df) <- NULL
  keggres.samedir.up.df <- as.data.frame(keggres.samedir$greater) ; keggres.samedir.up.df$dir <- "up" ; keggres.samedir.up.df$name <- rownames(keggres.samedir.up.df) ; rownames(keggres.samedir.up.df) <- NULL
  keggres.samedir.dn.df <- as.data.frame(keggres.samedir$less) ; keggres.samedir.dn.df$dir <- "down" ; keggres.samedir.dn.df$name <- rownames(keggres.samedir.dn.df) ; rownames(keggres.samedir.dn.df) <- NULL
  keggres.df <- rbind(keggres.df, keggres.samedir.up.df, keggres.samedir.dn.df)
  keggres.df <- keggres.df[order(keggres.df$q.val),]
  keggres.df <- keggres.df[!duplicated(keggres.df$name),] 
  keggres.df <- keggres.df[,c("name", colnames(keggres.df)[colnames(keggres.df)!="name"])]
  
  CSV.write(pathways.out, keggres.df, first.cell = "name")

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
  
  # draw pathway
  library(pathview)
  
  sig <- !is.na(keggres.df$q.val) & keggres.df$q.val <= fdrCutoff
  if (sum(sig) > 0) {
    for (pw in keggres.df$name[sig][1:min(length(keggres.df$name[sig]), topN)]) {
      enriched.name <- substr(pw, 1, 8)
      print(paste("Plotting pathway", enriched.name))
      enriched.id <- gsub("hsa", "", enriched.name)
      pw.file <- paste0(enriched.name, ".pathview.png")
      plot.file <- sprintf('%s-%s.png', instance.name, enriched.name)
      
      setwd(dirname(report.dir))
      pv.out <- pathview(gene.data = foldchanges, pathway.id = enriched.id, species = "hsa", low=list(gene="#6666FF"), mid=list(gene="#CCCCCC"), high=list(gene="#FF0000"), kegg.dir = dirname(report.dir))
      system(paste0("mv ", getwd(), "/", pw.file, " ", report.dir, "/", plot.file))
      if (!file.exists(paste0(report.dir, "/", plot.file))) {
        print(paste0("WARNING: Could not create pathway diagram for '", enriched.name, "'. Skipped!"))
        next
      }
      
      caption.prefix <- paste0("KEGG pathway '", pw, "'.")
      if (keggres.df$dir[keggres.df$name==pw] == "up") {
        caption.suffix <- sprintf("Highest enrichment (q=%.2g) was found in up-regulated genes.", keggres.df$q.val[keggres.df$name==pw])
      } else if (keggres.df$dir[keggres.df$name==pw] == "down") {
        caption.suffix <- sprintf("Highest enrichment (q=%.2g) was found in down-regulated genes.", keggres.df$q.val[keggres.df$name==pw])
      } else {
        caption.suffix <- sprintf("Highest enrichment (q=%.2g) was found considering both up- and down-regulated genes.", keggres.df$q.val[keggres.df$name==pw])
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

