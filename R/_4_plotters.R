## plot functions
message("import ploting module...")
qcPlots <- function(dataset, imgpath, addname = NULL, type, subset = "all", ref_levels = NA, cols="", col = "rev", use=c("display", "ppt", "publish", "poster"), title="") {
  
  # to generate qc plots.
  # including count density, cedf, sample distance heatmap, PCA, gene clusering heatmap

  ## check plot type and input object class
  if (!type %in% c("density", "ecdf", "heatmap", "PCA", "genecluster", "dispersion", "MA")) {
    stop("plot type unknown")
  }
  
  if (!is.na(ref_levels[1])) {
    dataset <- ddsFactorRelevel(dataset, ref_levels = ref_levels)
  }
  
  if (type %in% c("density", "ecdf")) {
    if (!class(dataset) == "DESeqDataSet") {
      stop("Wrong data type")
    }
    ## gene counts table
    geneCounts <- counts(dataset)
    keep <- rowSums(geneCounts) > ncol(geneCounts) * 10
    
  }
  
  if (type %in% c("heatmap", "PCA", "genecluster")) {
    if (class(dataset) == "DESeqDataSet") {
      cat("DESeqDataSet detected; performaing vst log2 transformation on dataset\n")
      cat("It may take a few minutes.\n")
      vst.dataset <- vst(dataset, blind = T)
    } else if (class(dataset) == "DESeqTransform") {
      cat("DESeqTransform datatype detected.\n")
    } else {
      stop("Wrong data type\n")
    }
  }
  
  #image export settings
  if (use == "ppt") {
    img <- list(type = "tiff", w = 10, h = 8, unit = "in", dpi = 150, px = 20, lsz = 1.5, psz = 6)
  } else if (use == "publish") {
    img <- list(type = "tiff", w = 5, h = 4, unit = "in", dpi = 300, px = 12)
  } else if (use == "poster" ) {
    img <- list(type = "tiff", w = 5, h = 4, unit = "in", dpi = 150)
  } else if (use == "display" ){
    img <- list(type = NA, w = 5, h = 4, unit = "in", dpi = 300, px = 12)
  } else {
    stop("unknown usage")
  }
  
  if (!use == "display") {
    #imgpath <- file.path(imgpath, "qcplots")
    prefix <- substring(rownames(colData(dataset))[1], 1,3)
    if (!"BCN" %in% levels(colData(dataset)$Treatment)) prefix <- paste0(prefix, "_CLE2")
    if (!is.null(addname)) prefix <- paste(prefix, addname, sep = "_")
    
    imgname <- paste(prefix, type, use, paste0(img$w, "x", img$h, img$unit), paste0(img$dpi,"dpi", paste0(".", img$type)), sep = "_")
    cat(imgpath)
    cat("\n")
    cat(imgname) 
  } else {
    imgpath <- NA
    imgname <- NA
  }

  fcts <- ddsGetFactors(dataset)
  dataset <- ddsCombineFactors(dataset, combined_name = "Group")
  
  fctlvls <- list()
  ncolors <- 1
  for (f in fcts) {
    fctlvls[[f]] <- levels(colData(dataset)[[f]])
    ncolors <- ncolors * length(fctlvls[[f]])
  }
  
  if(length(fctlvls[[1]])==2) {
    samplecol <- brewer.pal(10, "Paired")[1:min(ncolors,10)]
  } else if (length(fctlvls[[2]])==2) {
    samplecol <- c(brewer.pal(10, "Paired")[c(TRUE, FALSE)], brewer.pal(10, "Paired")[c(FALSE,TRUE)])[1:min(ncolors,10)]
  } else {
    samplecol <- brewer.pal(8, "Set2")[1:min(ncolors,8)]
  }
  if(col=="rev") {
    samplecol <- rev(samplecol)
  }
  labelcol <- colorMatch(colData(dataset)$Group, samplecol)
  
  ## set colors for plots
  hmcol <- colorRampPalette(brewer.pal(9, "OrRd"))(255)          ## color for heatmap
  gccol <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))  ## color for genecluster
  
  ## ploting
  if (type == "density" ) {
    #tiff(paste0(imgpath, imgname), width = img$w, height = img$h, units = img$unit, res = img$dpi, pointsize = img$px, compression = "lzw")
    multidensity(counts(dataset, normalized = T)[keep,], xlab="Normalized Counts", xlim=c(0,1000), col = labelcol, main = title)
    #dev.off()
    
  } else if (type == "ecdf" ) {
    multiecdf(counts(dataset, normalized = T)[keep,], xlab="Normalized Counts", xlim=c(0,1000), col = labelcol, main = title)
    
  } else if (type == "heatmap") {
    dataset <- ddsCombineFactors(dataset, combined_name = "Group")
    distsRL <- dist(t(assay(vst.dataset)))
    mat<- as.matrix(distsRL)
    colnames(mat) <- colData(dataset)$Group
    rownames(mat) <- rownames(colData(dataset))
    
    pheatmap(mat, trace="none", 
             fontsize = img$px, angle_col = 90, 
             col=rev(hmcol), border_color = NA, 
             filename = ifelse(is.na(imgname), NA, file.path(imgpath, imgname)), 
             width = img$w, height = img$w * 0.9,
             main=title)
    
  } else if (type == "PCA") {
    
    if (subset == "all") {
      ntop=nrow(dataset)
      tisufix <- "all genes"
    } else {
      ntop=subset
      tisufix <- sprintf("top variable genes (%s)", ntop)
    }
    
    ti <- sprintf("PC1 vs PC2, %s", tisufix)
    
    Pvars<- rowVars(assay(vst.dataset))
    select <- order(Pvars, decreasing=T)[seq_len(min(ntop, length(Pvars)))]
    PCA <- prcomp(t(assay(vst.dataset)[select,]), scale = F)
    percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
    
    dataGG <- data.frame(PC1=PCA$x[,1], PC2=PCA$x[,2], PC3=PCA$x[,3], PC4=PCA$x[,4],
                         sampleNo = rownames(colData(dataset)),
                         sample = colData(dataset)$Group)
    
    PCAplot <- ggplot(dataGG, aes(PC1, PC2, label = colnames(dataset))) +
              ggtitle(title) + 
              geom_point(aes(col=factor(sample)),size=5*img$px/14) +
              scale_colour_manual(values = samplecol) + 
              geom_text_repel(size = 3*img$px/14,
                              color=colRamp(labelcol,1.5),
                              hjust="inward", vjust = "inward", 
                              show.legend = F) +
              labs(title = ti,
                   x = paste0("PC1, VarExp:", round(percentVar[1], 4)),
                   y = paste0("PC2, VarExp:", round(percentVar[2], 4))) + 
              
              theme(#text = element_text(size = img$px),
                    axis.title.x = element_text(colour = "black",size = img$px),   # hide x axis title.
                    axis.ticks.x = element_line(size = img$lsz * 0.8, color = "black"), #hide x axis ticks.
                    axis.title.y = element_text(color = "black", size = img$px),
                    axis.ticks.y = element_line(size = img$lsz * 0.8, color = "black"),
                    axis.text.x = element_text(vjust = 1, hjust = 1, color = "black", size = img$px),
                    axis.text.y = element_text(face = "plain", color = "black", size = img$px),
                    
                    legend.background = element_blank(),
                    legend.key = element_rect(fill="white", color="white"),
                    legend.key.size = unit(img$px*0.8,"points"),
                    legend.text = element_text(size = img$px * 0.8, color = "black"),
                    legend.spacing.y = unit(img$px/5, "points"),
                    legend.title = element_blank(),
                    
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "white"),
                    panel.border = element_rect(color = "black", fill = NA, size = img$lsz),
                    plot.background = element_rect(fill = "white"),
                    plot.title = element_text(hjust = 0, size = img$px))

    if (!use=="display") {
        imgname <- sub(type, paste0(type, "(", ntop, ")"), imgname)
        ggsave(file.path(imgpath, imgname), device = img$type, scale = 1, width = img$w, height = img$w * 0.8, units = img$unit, dpi = img$dpi)   
      }
    
    return(PCAplot)
    
  } else if( type == "genecluster" ) {
        
        if (subset == "all") {
          ntop = nrow(dataset)
        } else {
          ntop = subset
        }
        
        topVar <- head(order(rowVars(counts(dataset, normalized=T)), decreasing = T),ntop)
        topVargenes <- as.matrix(counts(dataset, normalized=T)[topVar,])
        labels <- sub("(EWR|LCM)_", "", colnames(dataset))
        colnames(topVargenes) <- paste0(colData(dataset)$sample, "(", labels, ")") 
        title=sprintf("cluster of top %s variable genes", ntop)
        imgname <- sub(type, paste0(type, "(", ntop, ")"), imgname)
        #heatmap.2(topVargenes,scale="row", trace="none",dendrogram="column", col = gccol)
        pheatmap(topVargenes, 
                 scale="row", 
                 trace="none", 
                 col = gccol, 
                 border_color = NA,
                 angle_col = 90, 
                 fontsize = img$px,  
                 show_rownames = F,  
                 filename = ifelse(is.na(imgname), NA, file.path(imgpath, imgname)), 
                 width = img$w, height = img$w,
                 main=title)
      }
  
}

geneCluster <- function(dds, genelist, sortMethod = "groupVar", cluster_rows = T, ntop = NA, main = "", img.path = NA, img.name = NA) {
  
  if (is.na(ntop)) { 
    ntop <- length(genelist)
  } else {
    ntop <- min(length(genelist), ntop)
  }
  gccol <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))  ## color for genecluster
  
  dds <- ddsCombineFactors(dds)
  dds <- dds[match(genelist, row.names(dds)), ]
  row.names(dds) <- ifelse(mcols(dds)$Gene_name=="", row.names(dds), paste0(row.names(dds), " (", mcols(dds)$Gene_name, ")"))
  
  ## calculate group variances
  d <- counts(dds, normalized = T)
  novars <- which(rowVars(d)==0)
  d[novars, ] <- d[novars,] + abs(rnorm(ncol(d) * length(novars), 0, 0.00001))
  d.melt <- melt(d)
  d.melt <- merge(d.melt, as.data.frame(colData(dds)), by.x = "Var2", by.y = 0, all.x=T)
  d.melt <- select(d.melt, Var1, Group, value)
  d.melt$Group <- as.factor(d.melt$Group)
  m <- tapply(d.melt$value, list(d.melt$Var1, d.melt$Group), mean)
  #m <- m[order(row.names(m)), drop = F]
  
  if (sortMethod == "var") {
     keep <- head(order(rowVars(d), decreasing = T),ntop)
  } else if (sortMethod == "sum") {
     keep <- head(order(rowSums(d), decreasing = T), ntop)
  } else if (sortMethod == "mean") {
     keep <- head(order(rowMeans(d), decreasing = T), ntop)
  } else if (sortMethod == "groupVar") {
     keep <- head(order(rowVars(m), decreasing = T), ntop)
  } else {
     keep <- NA
  }
  
  if (!is.na(keep[1])) {
     mat <- as.matrix(d)[keep, , drop=F]
  } else {
     mat <- as.matrix(d)
  }

  if(nrow(mat) < 2 ) {
    message("Nothing to cluster: Less then 2 genes")
    return()
  }

  colnames(mat) <- dds$Group
  title <- main
  pheatmap(mat, 
           scale="row", 
           trace="none", 
           cluster_rows = cluster_rows,
           cellwidth = 60,
           cellheight = 12,
           col = gccol, 
           border_color = NA,
           angle_col = 90, 
           fontsize = 12,  
           show_rownames = T,  
           filename = ifelse(is.na(img.name), NA, file.path(img.path, img.name)), 
           #width = img$w, height = img$w,
           main = main)
  
}

plotCounts2 <- function(dds, genelist, intgroup="Genotype", facet = NULL, 
                        pty = "point", use.log = T, ncol = NA, cols = NULL,
                        countthreshold = 0, showTrend = T, showMean = T, 
                        showErrBar = F, legend = "right", font = 12,
                        main="", xlab="", ylab="Normalized Counts") {

  if (length(genelist) == 0) stop("no gene for ploting")
  if(ylab == "") ylab <- "Normalized Counts"
  mincount <- ncol(dds) * countthreshold
  keep <- which(rowSums(counts(dds, normalized = T)) >= mincount)
  dds <- dds[keep, ]
  genelist <- genelist[which(genelist %in% row.names(dds))]
  dds <- ddsCombineFactors(dds)
  fts <- ddsGetFactors(dds) ## extract factors from designs
  
  if(use.log) {
    # d.vst <- vst(dds, blind = T)
    # cnts <- assay(d.vst)
    cnts <- log2(counts(dds, normalized = T) + 1)
    ylab <- "log2(Expr)"
  } else {
    cnts <- counts(dds, normalized = T)
  }
  
  colnames(cnts) <- dds$Group
  tar_idx <- match(genelist, row.names(dds))
  if(all(is.na(tar_idx))) {
    stop("no gene to plot")
    return(0)
  }
  
  tar_idx <- tar_idx[!is.na(tar_idx)]
  cnts <- cnts[tar_idx, , drop = FALSE]
  #message(sprintf("count table dim:%s", str_c(dim(cnts), collapse = "\t")))
  #message(sprintf("count table rownames:%s", str_c(row.names(cnts), collapse = "\t")))
  gene_names_idx <- match(genelist, GENES_FOR_PLOT$Gene_id)
  gene_names <- GENES_FOR_PLOT$Gene_name[gene_names_idx]
  if(any(is.na(gene_names_idx))) {
    na_idx <- which(is.na(gene_names_idx))
    #gene_names[na_idx] <- mcols(dds)$Gene_name[match(genelist[na_idx], row.names(dds))]
    na_genes_idx <- match(genelist[na_idx], row.names(dds))
    gene_names[na_idx] <- as.character(mcols(dds)$Gene_name[na_genes_idx])
  }
  gene_names[gene_names == ""] <- NA
  
  ids <- row.names(cnts)
  row.names(cnts) <- sapply(1:nrow(cnts), 
                               function(i) {
                                  ifelse(is.na(gene_names[i]), 
                                          ids[i], 
                                          sprintf("%s (%s)", ids[i], gene_names[i]))
                                  })
  # row.names(cnts) <- ifelse(is.na(gene_names), 
                            # row.names(cnts), 
                            # sprintf("%s (%s)", row.names(cnts), gene_names))
  
  if(pty == "heatmap") {
    gccol <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(255))  ## color for genecluster
    plt <- pheatmap(cnts,
                     scale="row", trace="none", cluster_rows = T,
                     #cellwidth = 60, cellheight = 8,
                     col = gccol,  border_color = NA,
                     angle_col = 90,  fontsize = font, show_rownames = T,  
                     filename = NA, main = main)
    return(plt)
  }
  
  cnts.melt <- melt(cnts, value.name = "Counts")
  conds <- str_split(cnts.melt[,2], "\\.")
  conds_df <- data.frame(matrix(unlist(conds), nrow = length(conds), byrow = T))
  colnames(conds_df) <- fts
  conds_df <- data.frame(lapply(conds_df, 
                                function(x) factor(x, levels = unique(x))))
  cnts.melt <- cbind(cnts.melt, conds_df)
  colnames(cnts.melt)[1] <- "Gene"
  
  if (is.na(ncol)) {
    ncol <- ncolOptimize(length(genelist))
  }
  
  if(is.null(facet)) {
    x <- "Gene"
    facet <- fts[!fts %in% intgroup]
    ncol <- 1
  } else {
    x <- fts[!fts %in% intgroup]
  }
  
  grid.col <- "."
  grid.row <- facet
  if (pty == "bar") {
    showTrend <- F
    showMean <- F
    showErrBar <- T
  } 
  
  
  plt <- gg_Plot(cnts.melt, y = "Counts", x = x, group_by = intgroup, 
                 plt_type = pty, showTrend = showTrend, showMean = showMean, 
                 showErrBar = showErrBar, font = font, legend = legend, cols = cols,
                 grid.row = grid.row, grid.col = grid.col, ftcol = ncol, 
                 graph.title = main, xlab = xlab, ylab = ylab)
  
  return(plt)
}

gg_Plot <- function(df, y, x, group_by = NULL, plt_type = "point", cols = NULL,
                    showTrend = T, showMean = T, showErrBar = F, font = 12,
                    grid.row = ".", grid.col = ".", ftcol = 1, legend = "right", 
                    graph.title = "", xlab = "", ylab = "") {
  
  dodge.pts <- 0.25  # * (n-1)
  
  p <- ggplot(df, aes_string(x = x, y = y, group = group_by, col = group_by), na.rm = TRUE) 
  if(plt_type == "point") {
    geom_type <- geom_point(aes_string(col=group_by),
                            position = position_dodge(dodge.pts),
                            size=2)
  } else if (plt_type == "bar") {
    geom_type <- geom_bar(aes_string(fill=group_by),
                          color = "gray10",
                          position = position_dodge(0.9),
                          stat= "summary",
                          fun.y = "mean")
  }
  
  if (showErrBar) {
    err <- geom_errorbar(aes_string(group = group_by),
                         position = position_dodge(0.9), # make sure the number is the same as the one in geom_bar()
                         color = "gray10",
                         stat = "summary",
                         #fun.y = "mean_se",
                         width = 0.5)
  } else {
    err <- NULL
  }
  
  if(!is.null(cols)) {
    col.fill <- scale_fill_manual(values = cols)
    col.line <- scale_color_manual(values = cols)
  } else {
    col.fill <- NULL
    col.line <- NULL
  }
  
  if (showTrend) {
    trend <- geom_line(aes_string(group = group_by),
                   col = "gray85",
                   lwd = 0.75,
                   linetype = 1,
                   position = position_dodge(dodge.pts),
                   stat = "summary",
                   fun.y = "mean")
  } else {
    trend <- NULL
  }
  
  if (showMean) {
    m <- geom_point(aes_string(x = x, y = y, group = group_by),
                    pch = 3,
                    size = 3,
                    position = position_dodge(dodge.pts),
                    stat = "summary",
                    fun.y = "mean")
  } else {
    m <- NULL
  }

  ti <- ggtitle(graph.title)
  
  if(grid.row =="." & grid.col == ".") {
    ft <- NULL
  } else {
    ft <- facet_wrap(reformulate(grid.col, grid.row), ncol=ftcol, scales = 'free_y')
  }
  yx <- scale_y_continuous(name = ylab,  
                           limits = function(x) x * 1.2,
                           expand = c(0,0))
  yxlim <- expand_limits(y=0)
  xx <- scale_x_discrete(name = xlab, expand = c(0,0.25)) 
  
  th <- theme_bw(base_size = font) +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        # theme(axis.title.x = element_blank(),   # hide x axis title.
              # axis.ticks.x = element_line(size = font/20, color = "black"), #hide x axis ticks.
              # axis.title.y = element_text(face = "bold", color = "black", size = font * 1.2),
              # axis.ticks.y = element_line(size = font/20, color = "black"),
              # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", size = font),
              # axis.text.y = element_text(face = "plain", color = "black", size = font),
              # axis.line = element_line(color = "black", size = font/20, linetype = 1, lineend = "square"),
              # strip.text.x = element_text(face = "plain", colour = "black", size = font * 0.9, angle = 0),
              # panel.grid.major = element_blank(),
              # panel.grid.minor = element_blank(),
              # panel.background = element_rect(fill = "white"),
              # panel.border = element_rect(color = "black", fill = NA, size = font/20),
              # plot.background = element_rect(fill = "white"),
              # plot.title = element_text(hjust = 0.5, size = font * 1.5, face = "bold", colour = "black")) 
  lgd <- theme(legend.position = legend,
               legend.text = element_text(colour="black", size = font, face = "plain"),
               legend.title = element_blank())
  myplot <- p + ti + trend + geom_type + m + err + ft + xx + yx + yxlim + th + lgd + col.fill + col.line
  return(myplot)
}

colRamp <- function(col, amp){

      newcol <- col2rgb(col)
      if(amp <=1) {
         newcol <- newcol * amp # you can use of course other values than 2. Higher values the darker the output.
      } else if (amp > 1) {
         newcol <- newcol/amp + (amp-1)*255/amp
      }
      newcol <- rgb(t(newcol), maxColorValue=255)
      return(newcol)

}

ncolOptimize <- function(n) {
  n <- floor(n)
  if(n == 0) return(0)
  
  div <- seq_len(n)
  fts <- div[n %% div == 0]
  #f <- fts[1:ceiling(length(fts)/2)]
  
  nrow <- fts[ceiling(length(fts)/2)] # f[length(f)]
  ncol <- n/nrow
  
  if(ncol > 6 & ncol/nrow >= 2) {
    n <- n+1
    ncol <- ncolOptimize(n)
  }
  
  return(ncol)
}

colorMatch <- function(samples, cols) {
  nrep <- ceiling(length(samples)/length(cols))
  if(nrep > 1) cols <- rep(cols, nrep)
  
  if(!class(samples) == "factor") samples <- factor(samples, levels=unique(samples))
  if(!class(cols) == "factor") cols <- factor(cols, levels=unique(cols))
  
  c.df <- data.frame(col=cols, lvls = as.numeric(cols))
  s.df <- data.frame(sample = samples, lvls = as.numeric(samples))
  
  s.df$col <- levels(cols)[s.df$lvls]

  return(s.df$col)
  
}

paretoh <- function(x, label = NULL, col = NA, rampCol = F,  addline.v = NA, addline.h = NA, title=NULL, xlabel = NULL, ylabel=NULL, img.path = NA, img.name = NA) {
  
  # take a named numeric vector and return a barplot in vertical version
  
  textlength <- max(nchar(names(x)))
  op <- par()
  col <- ifelse(is.na(col), colRamp("red3",1.4), col)
  col.step <- x/max(x)
  if(rampCol)  col <- sapply(col.step, function(cs) colRamp(col, cs))
  
  if(!is.na(img.path)) {
      png(file.path(img.path, img.name), 
          width = textlength * 0.2, 
          height = length(x) * 0.3 + 1, 
          res = 300, unit = "in")
  }
  par(mar=c(4,textlength * 0.5,1,3))
  bp <- barplot(x, #names.arg = name, cex.names = 0.7,
                horiz = T,
                width = 0.2, space = 0.2, border = NA, axes = F,
                col = col,
                xlab = xlabel,
                cex.lab = 1.4,
                ylab = "",
                yaxt = "n",
                main = title)
  tt <- text(x = x, y = bp, label = label, pos = 4, cex = 0.75, col = "gray15", xpd = T)## add labels to top of bar
  ax2 <- axis(side = 2, at = bp, labels = names(x), tick = F, las= 2, xpd = F, col.axis = "gray15", col = "gray15", mgp=c(3,0,0), cex.axis=1.2)
  ax1 <- axis(side = 1, labels= T, col.axis = "gray15", col = "gray15", cex.axis = 1.2, las=1, mgp=c(1,1,0))
  if(!is.na(addline.h)) abline(h = addline.h, lty = 2, col = "gray90", lwd = 1.5)
  if(!is.na(addline.v)) abline(v = addline.v, lty = 2, col = "gray90", lwd = 1.5)
  abline(v=0, col = "gray15", lwd = 1.1)
  if (!is.na(img.path)) dev.off()
  suppressWarnings(par(op))
  
}

generateColors <- function(n, theme = "bw") {
     
    theme <- tolower(theme)
    if(theme %in% c("r", "red", "rd")) {
       theme <- "r"
    } else if(theme %in% c("blue", "bl", "b")) {
       theme <- "b"
    } else if(theme %in% c("green", "g")) {
       theme <- "g"
    } else if(theme %in% c("yellow", "y", "yl")) {
       theme <- "gr"
    } else if(theme %in% c("cyan", "c")) {
       theme <- "bg"
    } else if(theme %in% c("magenta", "m")) {
       theme <- "br"
    }
    rgb <- lapply(1:3, function(i) floor(seq(10,  225, length.out = n)))
    rgb <- lapply(rgb, as.hexmode)
    rgb <- lapply(rgb, as.character)
    rgb.df <- as.data.frame(rgb)
    colnames(rgb.df) <- c("r", "g", "b")   
    if(theme == "bw") {
       colors <- sprintf("#%s", apply(rgb.df, 1, str_c, collapse = ""))
    } else if (theme %in% c("r", "g", "b")) {
       rgb.df[, !colnames(rgb.df) == theme] <- "00"
       colors <- sprintf("#%s", apply(rgb.df, 1, str_c, collapse = ""))
    } else if (theme %in% c("gr", "br", "bg")) {
       FFcols <- sapply(colnames(rgb.df), function(x) grepl(x, theme))
       rgb.df[, FFcols] <- "FF"
       colors <- sprintf("#%s", apply(rgb.df, 1, str_c, collapse = ""))
    }
    
    return(colors)
    
}
