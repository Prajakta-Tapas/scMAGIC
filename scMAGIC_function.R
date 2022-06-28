library(Matrix)
library(data.table)
library(Seurat)
library(dplyr)
library(gt)
library(igraph)
library(RColorBrewer)
library(scMAGIC)

########################################
#Author: Yu Zhang
#Email: zhang_yu18@fudan.edu.cn
#######################################

get_overlap_genes <- function(exp_sc_mat, exp_ref_mat) {
  exp_ref_mat <- as.data.frame(exp_ref_mat)
  exp_sc_mat <- as.data.frame(exp_sc_mat)
  # get overlap genes
  exp_sc_mat <- exp_sc_mat[order(rownames(exp_sc_mat)),]
  exp_ref_mat <- exp_ref_mat[order(rownames(exp_ref_mat)),]
  gene_sc <- rownames(exp_sc_mat)
  gene_ref <- rownames(exp_ref_mat)
  gene_over <- gene_sc[which(gene_sc %in% gene_ref)]
  exp_sc_mat <- exp_sc_mat[gene_over,]
  exp_ref_mat <- exp_ref_mat[gene_over,]
  
  out.overlap <- list()
  out.overlap$exp_sc_mat <- exp_sc_mat
  out.overlap$exp_ref_mat <- exp_ref_mat
  out.overlap$gene_over <- gene_over
  return(out.overlap)
}


.get_high_variance_genes <- function(exp_ref_mat, num.genes = 2000, type_ref = 'sum-counts') {
  if (type_ref %in% c('sum-counts', 'sc-counts')) {
    seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
    seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                scale.factor = 1e6, verbose = F)
  }
  if (type_ref %in% c('fpkm', 'tpm', 'rpkm')) {
    exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
    seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
    seurat.Ref@assays$RNA@data <- exp_ref_mat
  }
  seurat.Ref <- FindVariableFeatures(
    seurat.Ref,
    selection.method = "vst",
    nfeatures = num.genes,
    verbose = F
  )
  return(VariableFeatures(seurat.Ref))
}


.one_multinomial <- function(i, exp_sc_mat, exp_ref_mat, colname_ref,
                             verbose, print_step) {
  delta <- 0.5
  Refprob <- function(exp_sc, exp_ref) {
    log_p_sc_given_ref <- dmultinom(x = exp_sc, log = T, prob = exp_ref)
    return(log_p_sc_given_ref)
  }
  #################
  exp_sc <- as.array(exp_sc_mat[, i])
  log_p_sc_given_ref_list <- numeric(length = length(colname_ref))
  j = 1
  while (j <= length(colname_ref)) {
    exp_ref <- as.array(exp_ref_mat[, j])
    #####
    exp_ref[which(exp_ref == 0)] <- delta * min(exp_ref[which(exp_ref > 0)])
    #####
    log_p_sc_given_ref <- Refprob(exp_sc, exp_ref)
    log_p_sc_given_ref_list[j] <- log_p_sc_given_ref
    j = j + 1
  }
  ################################
  if (verbose) {
    if (i %% print_step == 1) {
      print(i)
    }
  }
  
  return(log_p_sc_given_ref_list)
  
}


.get_log_p_sc_given_ref <- function(exp_sc_mat, exp_ref_mat, num_threads = 4,
                                    print_step = 100, gene_overlap = FALSE,
                                    verbose = FALSE) {
  library(parallel, verbose = F)
  
  #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
  #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name
  
  #Step 1. get overlapped genes
  if (gene_overlap) {
    if (verbose) {
      print('Gene number of exp_sc_mat:')
      print(nrow(exp_sc_mat))
      print('Gene number of exp_ref_mat:')
      print(nrow(exp_ref_mat))
    }
    out.overlap <- get_overlap_genes(exp_sc_mat, exp_ref_mat)
    exp_sc_mat <- out.overlap$exp_sc_mat
    exp_ref_mat <- out.overlap$exp_ref_mat
    if (verbose) {
      print('Number of overlapped genes:')
      print(nrow(exp_sc_mat))
    }
  }
  
  ###############
  colname_sc <- colnames(exp_sc_mat)
  colname_ref <- colnames(exp_ref_mat)
  
  #Step 2. calculate prob
  cl = makeCluster(num_threads, outfile = '')
  RUN <- parLapply(
    cl = cl,
    1:length(exp_sc_mat[1,]),
    .one_multinomial,
    exp_sc_mat = exp_sc_mat,
    exp_ref_mat = exp_ref_mat,
    colname_ref = colname_ref,
    verbose = verbose,
    print_step = print_step
  )
  stopCluster(cl)
  #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=num_threads)
  LOG_P_SC_GIVEN_REF = c()
  for (log_p_sc_given_ref_list in RUN) {
    LOG_P_SC_GIVEN_REF <-
      cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)
  }
  #######################################
  rownames(LOG_P_SC_GIVEN_REF) <- colname_ref
  colnames(LOG_P_SC_GIVEN_REF) <- colname_sc
  ######2019.02.16 start ######
  LOG_P_SC_GIVEN_REF[which(is.na(LOG_P_SC_GIVEN_REF))] <- min(LOG_P_SC_GIVEN_REF)
  ######2019.02.16 end ######
  return(LOG_P_SC_GIVEN_REF)
  
}


.one_get_corr <- function(barcode, exp_sc_mat, exp_ref_mat, colname_ref,
                          method, verbose, print_step) {
  # calculate a correlation
  exp_sc <- as.array(exp_sc_mat[, barcode])
  log_p_sc_given_ref_list <- numeric(length = length(colname_ref))
  j <- 1
  while (j <= length(colname_ref)) {
    exp_ref <- as.array(exp_ref_mat[, j])
    if (method == 'kendall') {
      log_p_sc_given_ref <- cor.fk(exp_sc, exp_ref)
    } else {
      if (method == 'cosine') {
        log_p_sc_given_ref <-
          sum(exp_sc * exp_ref) / sqrt(sum(exp_sc ^ 2) * sum(exp_ref ^ 2))
      } else {
        log_p_sc_given_ref <- cor(exp_sc, exp_ref, method = method)
      }
    }
    log_p_sc_given_ref_list[j] <- log_p_sc_given_ref
    j <- j + 1
  }
  ################################
  if (verbose) {
    if (i %% print_step == 1) {
      print(i)
    }
  }
  
  gc()
  return(log_p_sc_given_ref_list)
  
}


.get_cor <- function(exp_sc_mat, exp_ref_mat, method = 'kendall', num_threads = 4,
                     print_step = 100, gene_overlap = FALSE, verbose = FALSE){
  #method = "pearson", "kendall", "spearman"
  #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
  #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name
  
  #################
  library(parallel, verbose = F)
  #Step 1. get overlapped genes
  if (gene_overlap) {
    if (verbose) {
      print('Gene number of exp_sc_mat:')
      print(nrow(exp_sc_mat))
      print('Gene number of exp_ref_mat:')
      print(nrow(exp_ref_mat))
    }
    out.overlap <- get_overlap_genes(exp_sc_mat, exp_ref_mat)
    exp_sc_mat <- out.overlap$exp_sc_mat
    exp_ref_mat <- out.overlap$exp_ref_mat
    if (verbose) {
      print('Number of overlapped genes:')
      print(nrow(exp_sc_mat))
    }
  }
  
  ###############
  colname_sc <- colnames(exp_sc_mat)
  colname_ref <- colnames(exp_ref_mat)
  exp_sc_mat <- as.matrix(exp_sc_mat)
  exp_ref_mat <- as.matrix(exp_ref_mat)
  #######################################
  
  #Step 2. calculate corr
  cl <- makeCluster(num_threads, outfile = '')
  clusterEvalQ(cl, library(pcaPP))
  RUN <- parLapply(
    cl = cl,
    dimnames(exp_sc_mat)[[2]],
    .one_get_corr,
    exp_sc_mat = exp_sc_mat,
    exp_ref_mat = exp_ref_mat,
    colname_ref = colname_ref,
    method = method,
    verbose = verbose,
    print_step = print_step
  )
  stopCluster(cl)
  #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=num_threads)
  LOG_P_SC_GIVEN_REF <- c()
  for (log_p_sc_given_ref_list in RUN) {
    LOG_P_SC_GIVEN_REF <-
      cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)
  }
  #######################################
  rownames(LOG_P_SC_GIVEN_REF) <- colname_ref
  colnames(LOG_P_SC_GIVEN_REF) <- colname_sc
  ######2019.02.16 start ######
  LOG_P_SC_GIVEN_REF[which(is.na(LOG_P_SC_GIVEN_REF))] <- -1
  ######2019.02.16 end ######
  
  return(LOG_P_SC_GIVEN_REF)
}


.get_tag_max <- function(P_REF_GIVEN_SC) {
  RN <- rownames(P_REF_GIVEN_SC)
  CN <- colnames(P_REF_GIVEN_SC)
  TAG <- cbind(CN, rep('NA', length(CN)))
  i <- 1
  while (i <= length(CN)) {
    this_rn_index <- which(P_REF_GIVEN_SC[, i] == max(P_REF_GIVEN_SC[, i]))[1]
    TAG[i, 2] <- RN[this_rn_index]
    i <- i + 1
  }
  colnames(TAG) <- c('cell_id', 'tag')
  return(TAG)
}


.combine_tags <- function(df.tags1, df.tags2) {
  # concat reference pval and local pval
  pvalue1 <- df.tags1[, c('scRef.tag', 'pvalue')]
  names(pvalue1) <- c('scRef.tag.1', 'pvalue.1')
  pvalue2 <- df.tags2[, c('scRef.tag', 'pvalue')]
  names(pvalue2) <- c('scRef.tag.2', 'pvalue.2')
  pvalue <- merge(pvalue1, pvalue2, by = 'row.names')
  row.names(pvalue) <- pvalue$Row.names
  pvalue$Row.names <- NULL
  
  # select more confident tag
  mtx.tag <- as.matrix(pvalue[, c('scRef.tag.1', 'scRef.tag.2')])
  mtx.pval <- as.matrix(pvalue[, c('pvalue.1', 'pvalue.2')])
  mtx.rank <- apply(mtx.pval, 1, rank, ties.method = "first")
  tag.final <-
    apply(as.array(1:dim(mtx.tag)[1]), 1, function(i) {
      mtx.tag[i, mtx.rank[1, i]]
    })
  pval.final <-
    apply(as.array(1:dim(mtx.pval)[1]), 1, function(i) {
      mtx.pval[i, mtx.rank[1, i]]
    })
  tag.final <-
    data.frame(
      scRef.tag = tag.final,
      pvalue = pval.final,
      row.names = dimnames(mtx.tag)[[1]],
      stringsAsFactors = F
    )
  
  OUT <- list()
  OUT$pvalue <- pvalue
  OUT$tag.final <- tag.final
  return(OUT)
  
}


generate_ref <- function(exp_sc_mat, TAG, min_cell = 1, M = 'SUM',
                         refnames = FALSE ){
  M <- M
  # print(M)
  min_cell <- min_cell
  refnames <- refnames
  exp_sc_mat <- exp_sc_mat
  TAG <- TAG
  NewRef <- c()
  TAG[, 2] <- as.character(TAG[, 2])
  if (refnames == FALSE) {
    refnames <- names(table(TAG[, 2]))
  }
  else{
    refnames <- refnames
  }
  outnames <- c()
  for (one in refnames) {
    this_col <- which(TAG[, 2] == one)
    if (length(this_col) >= min_cell) {
      outnames <- c(outnames, one)
      if (length(this_col) > 1) {
        if (M == 'SUM') {
          this_new_ref <- apply(exp_sc_mat[, this_col], 1, sum)
        } else{
          this_new_ref <- apply(exp_sc_mat[, this_col], 1, mean)
        }
      }
      else{
        this_new_ref <- exp_sc_mat[, this_col]
      }
      NewRef <- cbind(NewRef, this_new_ref)
    }
  }
  if (is.null(dim(NewRef))) {
    return(NULL)
  }
  rownames(NewRef) <- rownames(exp_sc_mat)
  colnames(NewRef) <- outnames
  return(NewRef)
}


.imoprt_outgroup <- function(out.group = 'MCA', normalization = T, use_RUVseq = T) {
  if (class(out.group)[1] %in% c("data.frame", "matrix")) {
    df.out.group <- out.group
  } else {
    if (class(out.group)[1] == "character") {
      if (out.group %in% c("MCA", "HCL")) {
        if (out.group == "MCA") {
          data(df.MCA)
          df.out.group <- df.MCA
        } else {
          data(df.HCL)
          df.out.group <- df.HCL
        }
      } else {
        if (file.exists(out.group)) {
          file.out.group <- out.group
          df.out.group <-
            read.table(file.out.group, header = T, row.names = 1,
                       sep = '\t', check.names = F)
        }
      }
    } else {
      stop('Error: incorrect input of outgroup')
    }
  }
  if (normalization) {
    library(Seurat, verbose = F)
    seurat.out.group <-
      CreateSeuratObject(counts = df.out.group, project = "out.group",
                         min.cells = 1, min.features = 5000)
    seurat.out.group <-
      NormalizeData(seurat.out.group, normalization.method = "LogNormalize",
                    scale.factor = 1e6, verbose = F)
    if (use_RUVseq) {
      # get stably expression genes
      seurat.out.group <- FindVariableFeatures(
        seurat.out.group, selection.method = "mvp",
        mean.cutoff = c(3, Inf), dispersion.cutoff = c(-0.05, 0.05), verbose = F)
    }
    
    return(seurat.out.group)
    
  } else {
    return(df.out.group)
  }
  
}


.find_markers_cell <- function(cell.ref, list_near_cell, exp_ref_mat.cell, exp_ref_label,
                               mtx.combat, percent.high.exp,
                               mtx.combat.use, topN, i) {
  library(Seurat, verbose = F)
  cell <- cell.ref[i]
  # print(cell)
  vec.cell <- mtx.combat[, paste0('Ref.', cell)]
  # vec.cell.high <- vec.cell[vec.cell > quantile(vec.cell, percent.high.exp)]
  # vec.ref <- exp_ref_mat[, cell]
  # vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, percent.high.exp)]
  # high expression genes
  # genes.high <- intersect(names(vec.cell.high), names(vec.ref.high))
  # genes.high <- names(vec.ref.high)
  genes.high <- rownames(exp_ref_mat.cell)
  
  # function
  .getDEgeneF <- function(esetm = NULL, group = NULL, pair = FALSE,
                          block = NULL, p_adj = "fdr", fpkm = T) {
    # limma function
    if (is.null(esetm)) {
      cat(
        "esetm: gene expression matrix",
        "group: factor: \"c\"/\"d\"",
        "pair: TRUE/FALSE*",
        "block: e.g.1 2 2 1 if paired; blank if not",
        "p_adj: p.adjust, fdr* ",
        "fpkm: TRUE/FALSE*",
        sep = "\n"
      )
    } else{
      library(limma, verbose = F)
      if (pair) {
        design <- model.matrix( ~ block + group)
      } else{
        design <- model.matrix( ~ group)
      }
      suppressWarnings(fit <- lmFit(esetm, design))
      if (fpkm) {
        suppressWarnings(fit <- eBayes(fit, trend = T, robust = T))
      } else{
        suppressWarnings(fit <- eBayes(fit))
      }
      x <- topTable(fit, number = nrow(esetm), adjust.method = p_adj,
                    coef = "group2")
      x <- x[!is.na(row.names(x)), ]
      x <- x[!duplicated(row.names(x)), ]
      return(x)
    }
  }
  
  # diff in local reference and negative reference
  seurat.Ref.cell <- CreateSeuratObject(counts = exp_ref_mat.cell, project = "Ref")
  seurat.Ref.cell <- NormalizeData(seurat.Ref.cell, verbose = F)
  seurat.Ref.cell@meta.data$original.label <- exp_ref_label
  markers <- FindMarkers(seurat.Ref.cell, ident.1 = cell, group.by = 'original.label',
                         only.pos = T, features = genes.high, min.cells.group = 1,
                         min.pct = 0.1, min.diff.pct = 0,
                         logfc.threshold = 0.4, verbose = F)
  # markers$p_val_fdr <- p.adjust(markers$p_val, method = 'fdr')
  genes.ref.cell <- row.names(markers[(markers$p_val_adj < 0.05),])
  genes.ref <- genes.ref.cell
  
  # genes.ref = genes.high
  n.neighbor <- 3
  i <- 1
  for (cell_near in rev(list_near_cell[[cell]][1:n.neighbor])) {
    sub_seurat <- subset(seurat.Ref.cell, subset = original.label %in% c(cell, cell_near))
    sub_markers <- FindMarkers(sub_seurat, ident.1 = cell, group.by = 'original.label',
                               only.pos = T, features = genes.ref.cell, min.cells.group = 1,
                               min.pct = 0.1, min.diff.pct = 0,
                               logfc.threshold = 0.25, verbose = F)
    sub_markers <- sub_markers[order(sub_markers$p_val),]
    genes.ref <- intersect(genes.ref, row.names(sub_markers[(sub_markers$p_val < 0.05),]))
    coef <- 0.5*(1+n.neighbor-i)
    if (length(genes.ref) < (coef*topN)) {
      sub_markers <- sub_markers[order(sub_markers$p_val),]
      genes.ref <- row.names(sub_markers)[1:(coef*topN)]
    }
    i <- i + 1
  }
  
  if (length(genes.ref) > (2*topN)) {
    sub_markers <- sub_markers[order(sub_markers$p_val),]
    genes.ref <- row.names(sub_markers)[1:(2*topN)]
  }
  use.genes <- genes.ref
  
  # diff in Atlas
  if (!is.null(mtx.combat.use)) {
    if (length(use.genes) > 4*topN) {
      mtx.limma <- cbind(mtx.combat.use, vec.cell)
      bool.atlas.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
      res.limma.MCA <- .getDEgeneF(mtx.limma[use.genes, ], bool.atlas.cell)
      res.limma.MCA <- res.limma.MCA[res.limma.MCA$logFC > 0,]
      genes.diff <- row.names(res.limma.MCA[(res.limma.MCA$P.Value < 0.001),])
      if (length(genes.diff) < (topN)) {
        res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
        genes.diff <- row.names(res.limma.MCA)[1:(topN)]
      }
      if (length(genes.diff) > (3*topN)) {
        res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
        genes.diff <- row.names(res.limma.MCA)[1:(3*topN)]
      }
    } else {
      genes.diff <- use.genes
    }
  }else {
    genes.diff <- use.genes
  }
  
  return(genes.diff)
  
}

.find_markers <- function(exp_ref_mat, exp_ref_mat.cell, exp_ref_label, seurat.out.group,
                          type_ref = 'sum-counts', use_RUVseq = T,
                          base.topN = 50, percent.high.exp = 0.8, num_threads = 6) {
  library(parallel, verbose = F)
  library(Seurat, verbose = F)
  # check parameters
  if (!is.null(base.topN)) {
    base.topN <- base.topN
  } else {
    stop('Error in finding markers: provide incorrect parameters')
  }
  
  # transform count to fpm
  if (type_ref == 'sum-counts') {
    seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
    seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                scale.factor = 1e6, verbose = F)
    exp_ref_mat <- as.matrix(seurat.Ref@assays$RNA@data)
  }
  if (type_ref %in% c('fpkm', 'tpm', 'rpkm', 'bulk')) {
    exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
  }
  cell.ref <- dimnames(exp_ref_mat)[[2]]
  
  if (!is.null(seurat.out.group)) {
    ###### regard a outgroup (e.g. MCA/HCA) as reference of DEG
    seurat.MCA <- seurat.out.group
    fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
    
    # overlap genes
    fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
    out.overlap <- get_overlap_genes(fpm.MCA, exp_ref_mat)
    fpm.MCA <- as.matrix(out.overlap$exp_sc_mat)
    exp_ref_mat <- as.matrix(out.overlap$exp_ref_mat)
    if (use_RUVseq) {
      gene_overlap <- out.overlap$gene_over
      SEG.MCA <- VariableFeatures(seurat.MCA)
      gene.constant <- intersect(gene_overlap, SEG.MCA)
    }
    # print('Number of overlapped genes:')
    # print(nrow(exp_ref_mat))
    
    cell.MCA <- dimnames(fpm.MCA)[[2]]
    # use RUVseq to remove batch effect
    mtx.in <- cbind(fpm.MCA, exp_ref_mat)
    names.mix <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
    dimnames(mtx.in)[[2]] <- names.mix
    if (use_RUVseq) {
      library(RUVSeq, verbose = F)
      seqRUVg <- RUVg(as.matrix(mtx.in), gene.constant, k=1, isLog = T)
      mtx.combat <- seqRUVg$normalizedCounts
    } else {
      mtx.combat <- mtx.in
    }
    mtx.combat.use <- mtx.combat[, paste0('MCA.', cell.MCA)]
  } else {
    mtx.combat <- exp_ref_mat
    names.mix <- paste0('Ref.', cell.ref)
    dimnames(mtx.combat)[[2]] <- names.mix
    mtx.combat.use <- NULL
  }
  
  mat_cor <- cor(exp_ref_mat)
  # library(pcaPP)
  # cor.fk(exp_ref_mat)
  list_near_cell <- list()
  for (cell in cell.ref) {
    vec_cor <- mat_cor[cell,]
    list_near_cell[[cell]] <- names(sort(vec_cor, decreasing = T))[2:min(length(cell.ref), 4)]
  }
  
  topN <- base.topN
  cl = makeCluster(num_threads, outfile = '')
  # clusterExport(cl, '.getDEgeneF')
  RUN <- parLapply(
    cl = cl,
    1:length(cell.ref),
    .find_markers_cell,
    cell.ref = cell.ref, list_near_cell = list_near_cell,
    exp_ref_mat.cell = exp_ref_mat.cell, exp_ref_label = exp_ref_label,
    mtx.combat = mtx.combat, percent.high.exp = percent.high.exp,
    mtx.combat.use = mtx.combat.use, topN = topN
  )
  stopCluster(cl)
  names(RUN) <- cell.ref
  
  out <- list()
  out[['list.cell.genes']] <- RUN
  out[['list_near_cell']] <- list_near_cell
  return(out)
  
}


.find_markers_sc_cell <- function(cell.ref, list_near_cell, mtx.combat, LocalRef.sum, percent.high.exp,
                                  select.exp, vec.tag1, list.localNeg, mtx.combat.use, topN, i) {
  library(Seurat, verbose = F)
  cell <- cell.ref[i]
  # print(cell)
  vec.cell <- mtx.combat[, paste0('Ref.', cell)]
  # vec.cell.high <- vec.cell[vec.cell > quantile(vec.cell, percent.high.exp)]
  # vec.ref <- LocalRef.sum[, cell]
  # vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, percent.high.exp)]
  # high expression genes
  genes.high <- names(select.exp)
  
  # function
  .getDEgeneF <- function(esetm = NULL, group = NULL, pair = FALSE,
                          block = NULL, p_adj = "fdr", fpkm = T) {
    # limma function
    if (is.null(esetm)) {
      cat(
        "esetm: gene expression matrix",
        "group: factor: \"c\"/\"d\"",
        "pair: TRUE/FALSE*",
        "block: e.g.1 2 2 1 if paired; blank if not",
        "p_adj: p.adjust, fdr* ",
        "fpkm: TRUE/FALSE*",
        sep = "\n"
      )
    } else{
      library(limma, verbose = F)
      if (pair) {
        design <- model.matrix( ~ block + group)
      } else{
        design <- model.matrix( ~ group)
      }
      suppressWarnings(fit <- lmFit(esetm, design))
      if (fpkm) {
        suppressWarnings(fit <- eBayes(fit, trend = T, robust = T))
      } else{
        suppressWarnings(fit <- eBayes(fit))
      }
      x <- topTable(fit, number = nrow(esetm), adjust.method = p_adj,
                    coef = "group2")
      x <- x[!is.na(row.names(x)), ]
      x <- x[!duplicated(row.names(x)), ]
      return(x)
    }
  }
  
  # diff in local reference and negative reference
  seurat.Ref.cell <- CreateSeuratObject(counts = select.exp, project = "Ref")
  seurat.Ref.cell <- NormalizeData(seurat.Ref.cell, verbose = F)
  seurat.Ref.cell@meta.data$original.label <- vec.tag1
  markers <- FindMarkers(seurat.Ref.cell, ident.1 = cell, group.by = 'original.label',
                         only.pos = T, features = genes.high, min.cells.group = 1,
                         min.pct = 0.1, min.diff.pct = 0,
                         logfc.threshold = 0.4, verbose = F)
  # markers$p_val_fdr <- p.adjust(markers$p_val, method = 'fdr')
  genes.ref.cell <- row.names(markers[(markers$p_val_adj < 0.05),])
  genes.ref <- genes.ref.cell
  
  n.neighbor <- 3
  i <- 1
  for (cell_near in rev(list_near_cell[[cell]][1:n.neighbor])) {
    sub_seurat <- subset(seurat.Ref.cell, subset = original.label %in% c(cell, cell_near))
    sub_markers <- FindMarkers(sub_seurat, ident.1 = cell, group.by = 'original.label',
                               only.pos = T, features = genes.ref.cell, min.cells.group = 1,
                               min.pct = 0.1, min.diff.pct = 0,
                               logfc.threshold = 0.25, verbose = F)
    sub_markers <- sub_markers[order(sub_markers$p_val),]
    genes.ref <- intersect(genes.ref, row.names(sub_markers[(sub_markers$p_val < 0.05),]))
    coef <- 0.5*(1+n.neighbor-i)
    if (length(genes.ref) < (coef*topN)) {
      sub_markers <- sub_markers[order(sub_markers$p_val),]
      genes.ref <- row.names(sub_markers)[1:(coef*topN)]
    }
    i <- i + 1
  }
  
  if (length(genes.ref) > (2*topN)) {
    sub_markers <- sub_markers[order(sub_markers$p_val),]
    genes.ref <- row.names(sub_markers)[1:(2*topN)]
  }
  use.genes <- genes.ref
  
  # diff in Atlas
  if (!is.null(mtx.combat.use)) {
    if (length(use.genes) > 4*topN) {
      mtx.limma <- cbind(mtx.combat.use, vec.cell)
      bool.atlas.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
      res.limma.MCA <- .getDEgeneF(mtx.limma[use.genes, ], bool.atlas.cell)
      res.limma.MCA <- res.limma.MCA[res.limma.MCA$logFC > 0,]
      genes.diff <- row.names(res.limma.MCA[(res.limma.MCA$P.Value < 0.001),])
      if (length(genes.diff) < (topN)) {
        res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
        genes.diff <- row.names(res.limma.MCA)[1:min(topN, nrow(res.limma.MCA))]
      }
      if (length(genes.diff) > (2*topN)) {
        res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
        genes.diff <- row.names(res.limma.MCA)[1:(2*topN)]
      }
    } else {
      genes.diff <- use.genes
    }
  } else {
    genes.diff <- use.genes
  }
  
  return(genes.diff)
  
}


.find_markers_sc <- function(select.exp, vec.tag1, LocalRef,
                             seurat.out.group, list.localNeg,
                             use_RUVseq = T, num_threads = 6,
                             base.topN = 50, percent.high.exp = 0.80) {
  library(parallel, verbose = F)
  library(Seurat, verbose = F)
  # check parameters
  if (!is.null(base.topN)) {
    base.topN <- base.topN
  } else {
    stop('Error in finding markers: provide incorrect parameters')
  }
  
  # transform count to fpm
  seurat.Ref.sum <- CreateSeuratObject(counts = LocalRef, project = "Ref")
  seurat.Ref.sum <- NormalizeData(seurat.Ref.sum, verbose = F)
  LocalRef.sum <- as.matrix(seurat.Ref.sum@assays$RNA@data)
  cell.ref <- dimnames(LocalRef.sum)[[2]]
  
  if (!is.null(seurat.out.group)) {
    ###### regard a outgroup (e.g. MCA/HCL) as reference of DEG
    seurat.MCA <- seurat.out.group
    fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
    
    # overlap genes
    fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
    out.overlap <- get_overlap_genes(fpm.MCA, LocalRef.sum)
    fpm.MCA <- as.matrix(out.overlap$exp_sc_mat)
    LocalRef.sum <- as.matrix(out.overlap$exp_ref_mat)
    if (use_RUVseq) {
      gene_overlap <- out.overlap$gene_over
      SEG.MCA <- VariableFeatures(seurat.MCA)
      gene.constant <- intersect(gene_overlap, SEG.MCA)
    }
    # print('Number of overlapped genes:')
    # print(nrow(exp_ref_mat))
    
    cell.MCA <- dimnames(fpm.MCA)[[2]]
    cell.ref <- dimnames(LocalRef.sum)[[2]]
    # use RUVseq to remove batch effect
    mtx.in <- cbind(fpm.MCA, LocalRef.sum)
    names.mix <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
    dimnames(mtx.in)[[2]] <- names.mix
    if (use_RUVseq) {
      library(RUVSeq, verbose = F)
      seqRUVg <- RUVg(as.matrix(mtx.in), gene.constant, k=3, isLog = T)
      mtx.combat <- seqRUVg$normalizedCounts
    } else {
      mtx.combat <- mtx.in
    }
    mtx.combat <- mtx.in
    mtx.combat.use <- mtx.combat[, paste0('MCA.', cell.MCA)]
  } else {
    mtx.combat <- LocalRef.sum
    names.mix <- paste0('Ref.', cell.ref)
    dimnames(mtx.combat)[[2]] <- names.mix
    mtx.combat.use <- NULL
  }
  
  mat_cor <- cor(LocalRef.sum)
  # library(pcaPP)
  # cor.fk(exp_ref_mat)
  list_near_cell <- list()
  for (cell in cell.ref) {
    vec_cor <- mat_cor[cell,]
    list_near_cell[[cell]] <- names(sort(vec_cor, decreasing = T))[2:min(length(cell.ref), 4)]
  }
  
  topN <- base.topN
  cl = makeCluster(num_threads, outfile = '')
  # clusterExport(cl, '.getDEgeneF')
  RUN <- parLapply(
    cl = cl,
    1:length(cell.ref),
    .find_markers_sc_cell,
    cell.ref = cell.ref, list_near_cell = list_near_cell, mtx.combat = mtx.combat,
    LocalRef.sum = LocalRef.sum, percent.high.exp = percent.high.exp,
    select.exp = select.exp, vec.tag1 = vec.tag1,
    list.localNeg = list.localNeg, mtx.combat.use = mtx.combat.use, topN = topN
  )
  stopCluster(cl)
  names(RUN) <- cell.ref
  
  out <- list()
  out[['list.cell.genes']] <- RUN
  out[['exp_ref_mat']] <- LocalRef.sum
  return(out)
  
}


.cluster_sc <- function(exp_sc_mat, cluster_num_pc = 75, cluster_resolution = 3) {
  # cluster
  library(Seurat, verbose = F)
  # data preparing
  seurat.unlabeled <- CreateSeuratObject(counts = exp_sc_mat)
  seurat.unlabeled <-
    NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize",
                  scale.factor = 10000, verbose = F)
  # print(seurat.unlabeled@assays$RNA@data[1:6,1:6])
  seurat.unlabeled <-
    FindVariableFeatures(seurat.unlabeled, selection.method = "vst",
                         nfeatures = 2000, verbose = F)
  seurat.unlabeled <- ScaleData(seurat.unlabeled, verbose = F)
  
  # PCA
  seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = cluster_num_pc, verbose = F)
  
  # cluster
  seurat.unlabeled <-
    FindNeighbors(seurat.unlabeled, reduction = "pca", dims = 1:cluster_num_pc,
                  nn.eps = 0.5, verbose = F)
  seurat.unlabeled <-
    FindClusters(seurat.unlabeled, resolution = cluster_resolution,
                 n.start = 10, n.iter = 100, verbose = F)
  
  out.cluster <-
    data.frame(
      cluster.id = as.character(seurat.unlabeled@meta.data$seurat_clusters),
      # original.label = seurat.unlabeled@meta.data$original.label,
      row.names = dimnames(seurat.unlabeled@assays$RNA@counts)[[2]]
    )
  list_out <- list()
  list_out$out.cluster <- out.cluster
  list_out$seurat.query <- seurat.unlabeled
  return(list_out)
  
}


.cluster_increase_speed <- function(exp_sc_mat, df.cluster, combine_num_cell = 5, num_threads = 5) {
  library(parallel, verbose = F)
  sc.genes <- row.names(exp_sc_mat)
  df.cluster <- as.matrix(df.cluster)
  cluster.ids <- as.character(unique(df.cluster))
  combine_num_cell <- combine_num_cell
  
  .merge.one.cluster <- function(cluster.id) {
    library(Seurat)
    # merge cells in one cluster
    cell.ids <- names(df.cluster[df.cluster[, 'cluster.id'] == cluster.id,])
    sub.exp <- exp_sc_mat[, cell.ids]
    # print(dim(sub.exp))
    
    sub.seurat <- CreateSeuratObject(counts = sub.exp)
    sub.seurat <-
      NormalizeData(sub.seurat, normalization.method = "LogNormalize",
                    scale.factor = 10000, verbose = F)
    # print(sub.seurat@assays$RNA@data[1:6,1:6])
    sub.seurat <-
      FindVariableFeatures(sub.seurat, selection.method = "disp",
                           nfeatures = 1000, verbose = F)
    # print('1')
    sub.seurat <- ScaleData(sub.seurat,
                            features = VariableFeatures(sub.seurat),
                            verbose = F)
    # print(head(sub.seurat@assays$RNA@var.features))
    # print('2')
    
    # PCA
    num.cell <- dim(sub.exp)[2]
    if (num.cell < 5) {
      stop('Error: You should provide a smaller resolution!')
    }
    sub.seurat <- RunPCA(sub.seurat, npcs = min(10, (num.cell-2)), verbose = F)
    sub.pca <- sub.seurat@reductions$pca@cell.embeddings
    
    num.clusters <- ceiling(num.cell / combine_num_cell)
    if (num.clusters < 2) {
      sub.dict <- data.frame(
        cell.id = dimnames(sub.pca)[[1]],
        cluster.level1 = rep(cluster.id, num.cell),
        cluster.level2 = rep('1', num.cell)
      )
    } else {
      res.cluster <-
        kmeans(sub.pca, centers = num.clusters, nstart = 20, iter.max = 50)
      # names
      sub.dict <- data.frame(
        cell.id = names(res.cluster$cluster),
        cluster.level1 = rep(cluster.id, num.cell),
        cluster.level2 = res.cluster$cluster
      )
    }
    sub.dict$cluster.merge.id <-
      paste(sub.dict$cluster.level1, sub.dict$cluster.level2, sep = '-')
    row.names(sub.dict) <- sub.dict$cell.id
    
    # function
    .generate_ref <- function(exp_sc_mat, TAG, min_cell = 1, M = 'SUM',
                              refnames = FALSE ){
      M <- M
      # print(M)
      min_cell <- min_cell
      refnames <- refnames
      exp_sc_mat <- exp_sc_mat
      TAG <- TAG
      NewRef <- c()
      TAG[, 2] <- as.character(TAG[, 2])
      if (refnames == FALSE) {
        refnames <- names(table(TAG[, 2]))
      }
      else{
        refnames <- refnames
      }
      outnames <- c()
      for (one in refnames) {
        this_col <- which(TAG[, 2] == one)
        if (length(this_col) >= min_cell) {
          outnames <- c(outnames, one)
          if (length(this_col) > 1) {
            if (M == 'SUM') {
              this_new_ref <- apply(exp_sc_mat[, this_col], 1, sum)
            } else{
              this_new_ref <- apply(exp_sc_mat[, this_col], 1, mean)
            }
          }
          else{
            this_new_ref <- exp_sc_mat[, this_col]
          }
          NewRef <- cbind(NewRef, this_new_ref)
        }
      }
      if (is.null(dim(NewRef))) {
        return(NULL)
      }
      rownames(NewRef) <- rownames(exp_sc_mat)
      colnames(NewRef) <- outnames
      return(NewRef)
    }
    
    # merge expression profile
    tag.in <- sub.dict[, c('cell.id', 'cluster.merge.id')]
    sub.exp.merge <- .generate_ref(sub.exp, tag.in)
    sub.exp.merge <- sub.exp.merge[sc.genes, ]
    # print(class(sub.exp.merge))
    if (class(sub.exp.merge)[1] %in% c('numeric', 'integer')) {
      sub.exp.merge <- as.data.frame(sub.exp.merge)
      names(sub.exp.merge) <- unique(sub.dict$cluster.merge.id)
    }
    
    sub.out <- list()
    sub.out$sub.dict <- sub.dict
    sub.out$sub.exp.merge <- sub.exp.merge
    return(sub.out)
    
  }
  
  
  # split dataset
  cl.input <- list()
  cl <- makeCluster(num_threads, outfile = '')
  # cl <- makeCluster(num_threads, outfile = '', type = 'FORK')
  # clusterExport(cl, '.generate_ref')
  # clusterEvalQ(cl, library(Seurat))
  out.par <- parLapply(
    cl = cl,
    cluster.ids,
    .merge.one.cluster
  )
  stopCluster(cl)
  
  df.dict <- data.frame(stringsAsFactors = F)
  df.exp.merge <- data.frame(stringsAsFactors = F)
  i = 1
  for (sub.out in out.par) {
    df.dict <- rbind(df.dict, sub.out$sub.dict)
    if (i == 1) {
      df.exp.merge <- sub.out$sub.exp.merge
    } else {
      df.exp.merge <- cbind(df.exp.merge, sub.out$sub.exp.merge)
    }
    i = i + 1
  }
  seurat.exp.merge <-
    CreateSeuratObject(counts = df.exp.merge, project = "exp.merge",
                       min.cells = 1, min.features = 1)
  seurat.exp.merge <-
    NormalizeData(seurat.exp.merge, normalization.method = "LogNormalize", verbose = F)
  df.exp.merge <- seurat.exp.merge@assays$RNA@data
  
  out.merge <- list()
  out.merge$df.dict <- df.dict
  out.merge$df.exp.merge <- df.exp.merge
  return(out.merge)
  
}


.confirm_label_auc <- function(exp_sc_mat, list.cell.genes, list_near_cell, scRef.tag, num_threads = 10) {
  library(parallel, verbose = F)
  library(AUCell, verbose = F)
  set.seed(123)
  # confirm label
  exp_sc_mat <- as.matrix(exp_sc_mat)
  df.tag <- as.data.frame(scRef.tag)
  row.names(df.tag) <- df.tag$cell_id
  df.tag$cell_id <- NULL
  names(df.tag) <- 'scRef.tag'
  
  total_genes <- unique(unlist(list.cell.genes))
  df.auc <- data.frame()
  list.auc.back <- list()
  for (cell in names(list.cell.genes)) {
    if (sum(df.tag$scRef.tag==cell) < 1) {
      next
    }
    sel_cols <- colnames(exp_sc_mat)[(colnames(exp_sc_mat)%in%(colnames(exp_sc_mat)[df.tag$scRef.tag==cell]))]
    total_mtx <- exp_sc_mat[total_genes, sel_cols]
    if (is.null(dim(total_mtx))) {
      total_mtx <- matrix(total_mtx, ncol=1)
      dimnames(total_mtx)[[1]] <- total_genes
      dimnames(total_mtx)[[2]] <- sel_cols
    } else {
      total_mtx <- as.matrix(total_mtx)
    }
    cells_rankings <- AUCell_buildRankings(total_mtx, verbose = F, plotStats = F)
    
    cell_markers <- list.cell.genes[[cell]]
    cells_near <- list_near_cell[[cell]]
    genes_near <- unique(c(list.cell.genes[[cells_near[1]]], list.cell.genes[[cells_near[2]]],
                           list.cell.genes[[cells_near[3]]]))
    # genes_near <- list.cell.genes[[cells_near[1]]]
    # genes_near <- unique(c(list.cell.genes[[cells_near[1]]], list.cell.genes[[cells_near[2]]]))
    num_genes <- min(length(cell_markers), length(genes_near))
    genes.marker <- sample(cell_markers, num_genes)
    genes.back <- sample(genes_near, num_genes)
    # genes.marker <- cell_markers
    # genes.back <- genes_near
    
    geneSets <- list(geneSet1=genes.marker)
    cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank = length(genes.marker)/2,
                                verbose = F)
    
    # df.tags1_cd8 <- df.tags1[colnames(exp_sc_mat)%in%
    #                              (colnames(exp_sc_mat)[tag1[,2]==cell]),]
    # df.tags1_cd8$auc <- unlist(cells_AUC@assays@data@listData)
    # View(df.tags1_cd8)
    # background
    geneSets.back <- list(geneSet1=genes.back)
    cells_AUC.back <- AUCell_calcAUC(geneSets.back, cells_rankings,
                                     aucMaxRank = length(genes.back)/2,
                                     verbose = F)
    list.auc.back[[cell]] <- unlist(cells_AUC.back@assays@data@listData)
    df.auc <- rbind(df.auc, data.frame(AUC = unlist(cells_AUC@assays@data@listData),
                                       AUC_back = unlist(cells_AUC.back@assays@data@listData),
                                       row.names = colnames(total_mtx)))
    
  }
  meta.tag <- merge(df.tag, df.auc, by = 'row.names')
  row.names(meta.tag) <- meta.tag$Row.names
  meta.tag$Row.names <- NULL
  # meta.tag$log10Pval <- -log10(meta.tag$pvalue)
  
  gc()
  
  list_out <- list()
  list_out$meta.tag <- meta.tag
  list_out$list.auc.back <- list.auc.back
  return(list_out)
  
}


.one_cutoff_AUC <- function(i, cells, df.tags1, auc_gap, list_tags1_back) {
  library(mclust, verbose = F)
  cell <- cells[i]
  df.sub <- df.tags1[df.tags1$scRef.tag == cell, ]
  sub_AUC <- df.sub[, c('AUC', 'diff')]
  cell_back <- list_tags1_back[[cell]]
  back_value <- min(quantile(cell_back,0.9), 0.05)
  if (sum(df.tags1$scRef.tag == cell) < 5) {
    one_out <- list()
    one_out[[1]] <- cell
    one_out[[2]] <- c(max(max(sub_AUC$AUC), back_value + 0.2),
                      max(max(sub_AUC$diff), 0.05))
    one_out[[3]] <- median(cell_back)
  } else {
    if (max(sub_AUC[,1]) - min(sub_AUC[,1]) > 0.57) {
      min_G <- 4
      max_G <- 6
    } else {
      min_G <- 2
      max_G <- 5
    }
    model_AUC <- densityMclust(sub_AUC[,1], verbose = F, G = min_G:max_G)
    mean_AUC <- model_AUC$parameters$mean
    if (length(table(model_AUC$classification)) != length(mean_AUC)) {
      mean_AUC <- mean_AUC[names(table(model_AUC$classification))]
    }
    diff_AUC <- c()
    for (i in 1:(length(mean_AUC)-1)) {
      diff_AUC <- c(diff_AUC, mean_AUC[i+1] - mean_AUC[i])
    }
    if (nrow(df.sub) < 200) {
      auc_gap <- max(auc_gap, 0.2)
    }
    if (sum(diff_AUC > auc_gap) > 0) {
      cut_class <- max(as.numeric(names(mean_AUC)[1:length(diff_AUC)])[diff_AUC > auc_gap])
      # cut_AUC <- mean_AUC[cut_class]
      df.sub$class_AUC <- model_AUC$classification
      cut_AUC <- max(df.sub$AUC[df.sub$class_AUC == as.character(cut_class)])
      # cut_AUC <- quantile(df.sub$AUC[df.sub$class_AUC == as.character(cut_class)], 0.95)
    } else {
      cut_AUC <- 0
    }
    
    sub_AUC_scale <- scale(sub_AUC)
    model <- densityMclust(sub_AUC_scale, verbose = F, G = min_G:max_G)
    # cutoff_pos <- max(summary(sub_AUC)[5], max(cell_back))
    # cutoff_pos <- max(quantile(sub_AUC, 0.85), quantile(cell_back,0.95))
    # cutoff_pos <- max(0.85*(max(sub_AUC)-min(sub_AUC))+min(sub_AUC), quantile(cell_back,0.95))
    df.sub$classification <- model$classification
    cluster.mean <- as.data.frame(model$parameters$mean)[, as.numeric(names(table(model$classification)))]
    sel_class <- names(table(model$classification))[order(colSums(cluster.mean), decreasing = T)[1]]
    if (sum(df.sub$classification == sel_class) <= 6) {
      first_two <- order(colSums(cluster.mean), decreasing = T)[1:2]
      first_two_diff <- mean(df.sub[df.sub$classification == first_two[1], 'AUC']) -
        mean(df.sub[df.sub$classification == first_two[2], 'AUC'])
      if (first_two_diff < 0.2) {
        sel_class <- names(table(model$classification))[order(colSums(cluster.mean), decreasing = T)[1:2]]
      }
    }
    df.sub_sel <- df.sub[df.sub$classification %in% sel_class,]
    one_out <- list()
    one_out[[1]] <- cell
    one_out[[2]] <- c(max(quantile(df.sub_sel$AUC, 0.25), back_value),
                      max(quantile(df.sub_sel$diff, 0.25), 0.05))
    one_out[[3]] <- cut_AUC
    # one_out[[3]] <- median(cell_back)
  }
  
  return(one_out)
}


.cutoff_AUC <- function(df.tags1, list_tags1_back, exp_sc_mat, threshold, num_threads = num_threads) {
  library(parallel, verbose = F)
  cells <- as.character(unique(df.tags1$scRef.tag))
  list.cutoff <- list()
  vec.neg.cutoff <- c()
  vec.cut <- c()
  if (median(colSums(exp_sc_mat != 0)) < 1200) {
    base_thre <- 0.23
  } else {
    base_thre <- 0.18
  }
  if (threshold <= 5) {
    auc_gap <- (5-threshold)*2 + base_thre
  } else {
    auc_gap <- (5-threshold) + base_thre
  }
  
  cl = makeCluster(num_threads, outfile = '')
  RUN <- parLapply(
    cl = cl,
    1:length(cells),
    .one_cutoff_AUC,
    cells = cells,
    df.tags1 = df.tags1,
    auc_gap = auc_gap,
    list_tags1_back = list_tags1_back
  )
  stopCluster(cl)
  for (one_out in RUN) {
    cell <- one_out[[1]]
    list.cutoff[[cell]] <- one_out[[2]]
    # list.cutoff <- c(quantile(df.sub_sel$AUC, 0.25), quantile(df.sub_sel$diff, 0.25))
    # vec.neg.cutoff <- c(vec.neg.cutoff, one_out[[3]])
    vec.cut <- c(vec.cut, one_out[[3]])
  }
  # names(vec.cutoff) <- cells
  # names(vec.neg.cutoff) <- cells
  names(vec.cut) <- cells
  
  out.cutoff <- list()
  out.cutoff$list.cutoff <- list.cutoff
  out.cutoff$vec.neg.cutoff <- vec.neg.cutoff
  out.cutoff$vec.cut <- vec.cut
  return(out.cutoff)
}


scMAGIC <- function(exp_sc_mat, exp_ref_mat, exp_ref_label = NULL,
                    single_round = F, identify_unassigned = T,
                    atlas = NULL, use_RUVseq = T,
                    cluster_num_pc = 50, cluster_resolution = 3,
                    combine_num_cell = NULL, min_cell = 1,
                    method1 = 'kendall', method2 = NULL,
                    corr_use_HVGene1 = 2000, corr_use_HVGene2 = 2000,
                    threshold = 5, num_threads = 4, cluster_assign = F,
                    simple_output = T) {
  
  library(parallel, verbose = F)
  library(scibet)
  # check parameters
  # if (!type_ref %in% c('sc-counts', 'sum-counts', 'fpkm', 'tpm', 'rpkm')) {
  #     stop('Error: inexistent input of reference data format')
  # }
  cutoff.1 = 'default'
  cutoff.2 = 'default'
  mod = ''
  # simple_output = T
  type_ref = 'sc-counts'
  out.group = atlas
  opt_speed = F
  num_cell <- ncol(exp_sc_mat)
  if (is.null(combine_num_cell)) {
    if (num_cell > 3000) {
      combine_num_cell = 5
    } else {
      combine_num_cell = 3
    }
  }
  if (is.null(method2)) {
    if (num_cell > 5000) {
      method2 = 'randomforest'
    } else {
      method2 = 'multinomial'
    }
  }
  
  time1 <- Sys.time()
  # get sum-counts format
  if (type_ref == 'sc-counts') {
    print('Sum single cell counts matrix:')
    label.in <- data.frame(cell_id = colnames(exp_ref_mat),
                           tag = as.character(exp_ref_label))
    exp_ref_mat.cell <- exp_ref_mat
    exp_ref_mat.sum <- generate_ref(exp_ref_mat, label.in, M = 'SUM')
    exp_ref_mat <- exp_ref_mat.sum
    type_ref <- 'sum-counts'
  }
  
  # get overlap genes
  out.overlap <- get_overlap_genes(exp_sc_mat, exp_ref_mat)
  exp_sc_mat <- out.overlap$exp_sc_mat
  exp_ref_mat <- out.overlap$exp_ref_mat
  gene_over <- out.overlap$gene_over
  # print('Number of overlapped genes:')
  # print(nrow(exp_sc_mat))
  
  if (!is.null(out.group)) {
    seurat.out.group <- .imoprt_outgroup(out.group = out.group, use_RUVseq = use_RUVseq)
    # overlap genes
    gene.overlap <- intersect(gene_over, rownames(seurat.out.group@assays$RNA@counts))
    exp_sc_mat <- exp_sc_mat[gene.overlap, ]
    exp_ref_mat <- exp_ref_mat[gene.overlap, ]
    exp_ref_mat.cell <- exp_ref_mat.cell[gene.overlap,]
  } else {
    seurat.out.group <- NULL
  }
  
  print('Number of overlapped genes:')
  print(nrow(exp_sc_mat))
  
  # cluster analysis
  print('Start clustering :')
  list_out.cluster <-
    .cluster_sc(exp_sc_mat,
                cluster_num_pc = cluster_num_pc,
                cluster_resolution = cluster_resolution)
  df.cluster <- list_out.cluster$out.cluster
  seurat.query <- list_out.cluster$seurat.query
  
  out.merge <-
    .cluster_increase_speed(exp_sc_mat, df.cluster,
                            combine_num_cell = combine_num_cell, num_threads = num_threads)
  df.dict <- out.merge$df.dict
  print('Clustering completed!')
  # speed calculation
  if (opt_speed) {
    print('Speed calculation by clustering:')
    # out.merge <-
    #     .cluster_increase_speed(exp_sc_mat, df.cluster,
    #                             combine_num_cell = combine_num_cell, num_threads = num_threads)
    # df.dict <- out.merge$df.dict
    df.exp.merge <- out.merge$df.exp.merge
    df.dict.merge <- unique(df.dict[,c('cluster.level1', 'cluster.merge.id')])
    rownames(df.dict.merge) <- df.dict.merge$cluster.merge.id
    df.dict.merge <- df.dict.merge[colnames(df.exp.merge),]
    min_cell <- ceiling(min_cell / combine_num_cell)
    query_set <- as.data.frame(t(df.exp.merge))/1.0
    query_set$label <- df.dict.merge$cluster.level1
  } else {
    df.exp.merge <- exp_sc_mat
    query_set <- as.data.frame(t(df.exp.merge))/1.0
    query_set$label <- df.cluster$cluster.id
  }
  
  gene_not0 <- rownames(df.exp.merge)[rowSums(df.exp.merge)!=0]
  df.exp.merge <- df.exp.merge[gene_not0, ]
  exp_ref_mat <- exp_ref_mat[gene_not0, ]
  exp_ref_mat.cell <- exp_ref_mat.cell[gene_not0, ]
  
  # find markers of cell types in reference
  topN = 50
  percent.high.exp = 0.8
  print('Find marker genes of cell types in reference:')
  suppressMessages(
    out.markers <-
      .find_markers(
        exp_ref_mat, exp_ref_mat.cell, exp_ref_label,
        seurat.out.group,
        type_ref = 'sum-counts',
        use_RUVseq = use_RUVseq,
        percent.high.exp = percent.high.exp, num_threads = num_threads
      ))
  list.cell.genes <- out.markers[['list.cell.genes']]
  list_near_cell <- out.markers[['list_near_cell']]
  df.exp.merge <- exp_sc_mat
  
  # rm(exp_sc_mat)
  gc()
  df.exp.merge <- as.matrix(df.exp.merge)
  
  print('First-round annotation:')
  # print(method1)
  if (method1 == 'xxx (based PC)') {
    # pred_tags <- func_sciBet(df.exp.merge, exp_ref_mat.cell, exp_ref_label)
    # tag1 <- data.frame(cell_id = colnames(df.exp.merge), tag = pred_tags)
    # tag1 <- as.matrix(tag1)
    
    # pred_tags <- func_scClassify(df.exp.merge, exp_ref_mat.cell, exp_ref_label)
    # # pred_tags=func_scPred(df.exp.merge, exp_ref_mat.cell, exp_ref_label)
    # # pred_tags=func_seurat(df.exp.merge, exp_ref_mat.cell, exp_ref_label)
    # table(label_sc, pred_tags)
  } else {
    if (!is.null(corr_use_HVGene1)) {
      num_feature <- corr_use_HVGene1/2
      train_set <- as.data.frame(t(exp_ref_mat.cell))/1.0
      train_set$label <- exp_ref_label
      HVG <- SelectGene(train_set, k = num_feature)
      similarity.in <- df.exp.merge[HVG, ]
      ref.in <- exp_ref_mat[HVG, ]
    } else {
      similarity.in <- df.exp.merge
      ref.in <- exp_ref_mat
    }
    if (method1 != 'multinomial') {
      out1 <- .get_cor(ref.in, similarity.in, method = method1, num_threads = num_threads)
      out1 <- t(out1)
      tag1 <- .get_tag_max(out1)
    } else {
      test_set <- query_set[, etest_gene]
      pred_tags <- SciBet(train_set, test_set, k=num_feature, result = 'list')
      out1 <- SciBet(train_set, test_set, k=num_feature, result = 'table')
      rownames(out1) <- rownames(test_set)
      out1 <- t(out1)
      tag1 <- data.frame(cell_id = rownames(test_set), tag = pred_tags)
      tag1 <- as.matrix(tag1)
    }
  }
  out1_scale <- scale(out1)
  out1_diff <- apply(out1_scale, 2, function(input) {sort(input, decreasing = T)[1] - sort(input, decreasing = T)[2]})
  
  gc()
  
  print('Build local reference')
  list_out_1 <- .confirm_label_auc(df.exp.merge, list.cell.genes, list_near_cell,
                                   tag1, num_threads = num_threads)
  df.tags1 <- list_out_1$meta.tag
  list_tags1_back <- list_out_1$list.auc.back
  cell_ids <- colnames(df.exp.merge)
  df.tags1 <- df.tags1[cell_ids, ]
  
  df.tags1$diff <- out1_diff
  df.tags1 <- merge(df.tags1, df.dict, by = 'row.names')
  rownames(df.tags1) <- df.tags1$Row.names
  df.tags1$Row.names <- NULL
  out.cutoff <- .cutoff_AUC(df.tags1, list_tags1_back, exp_sc_mat, threshold, num_threads = num_threads)
  df.cutoff.1 <- out.cutoff$list.cutoff
  neg.cutoff.1 <- out.cutoff$vec.neg.cutoff
  vec.cut_1 <- out.cutoff$vec.cut
  
  if (single_round) {
    pvalue1 <- data.frame(stringsAsFactors = F)
    for (one_cell in unique(df.tags1$scRef.tag)) {
      df_sub <- df.tags1[df.tags1$scRef.tag == one_cell, c('scRef.tag', 'AUC')]
      df_sub$tag_unassaigned <- df_sub$scRef.tag
      df_sub$tag_unassaigned[df_sub$AUC < vec.cut_1[one_cell]] <- 'Unassigned'
      pvalue1 <- rbind(pvalue1, df_sub)
    }
    names(pvalue1) <- c('scRef.tag.1', 'AUC.1', 'tag_unassaigned.1')
    df_unassigned <- pvalue1[cell_ids, ]
    df_unassigned$tag.final <- pvalue1$scRef.tag.1
    df_unassigned$tag.final[df_unassigned$tag_unassaigned.1 == 'Unassigned'] <- 'Unassigned'
    
    df.tags <- merge(df_unassigned, df.dict, by = 'row.names')
    row.names(df.tags) <- df.tags$Row.names
    df.tags$Row.names <- NULL
    df.tags <- df.tags[, c("tag.final", "cluster.merge.id", "cluster.level1")]
    df.tags$scRef.tag <- df.tags$tag.final
    all_sub_clusters <- unique(df.tags$cluster.merge.id)
    for (sub_cluster in all_sub_clusters) {
      df_sub <- df.tags[df.tags$cluster.merge.id == sub_cluster,]
      sub_table <- sort(table(df_sub$tag.final), decreasing = T)
      if (length(sub_table) == 1) {
        next
      } else {
        if (sub_table[1]/nrow(df_sub) >= 0.6) {
          df.tags$scRef.tag[df.tags$cluster.merge.id == sub_cluster] <- names(sub_table)[1]
        }
      }
    }
    
    cell_ids <- colnames(exp_sc_mat)
    df.tags <- df.tags[cell_ids, ]
    
    if (cluster_assign) {
      all_clusters <- unique(df.tags$cluster.level1)
      df.tags$cluster.tags <- df.tags$scRef.tag
      for (cluster in all_clusters) {
        df_sub <- df.tags[df.tags$cluster.level1 == cluster,]
        sub_table <- sort(table(df_sub$tag.final), decreasing = T)
        df.tags$cluster.tags[df.tags$cluster.level1 == cluster] <- names(sub_table)[1]
      }
      df.combine <- data.frame(scMAGIC.tag = df.tags$cluster.tags, row.names = rownames(df.tags))
    } else {
      df.combine <- data.frame(scMAGIC.tag = df.tags$scRef.tag, row.names = rownames(df.tags))
    }
    
    gc()
    
    time2 <- Sys.time()
    time.scRef <- difftime(time2, time1, units = 'secs')
    output <- list()
    output$tag1 <- tag1
    output$out1 <- out1
    output$combine.out <- df.tags
    output$dict.cluster <- df.dict
    output$ref.markers <- list.cell.genes
    output$final.out <- df.combine
    output$run.time <- time.scRef
    
    print('Finish!')
    
    return(output)
    
  }
  
  select.barcode <- c()
  for (cell in names(df.cutoff.1)) {
    sub.cutoff <- df.cutoff.1[[cell]]
    sub.select <- df.tags1[df.tags1$scRef.tag == cell, ]
    sub.select <- sub.select[sub.select$AUC >= sub.cutoff[1] & sub.select$diff >= sub.cutoff[2], ]
    select.barcode <- c(select.barcode, row.names(sub.select))
  }
  list.localNeg <- list()
  
  df.tags1$conf <- rep(0, nrow(df.tags1))
  df.tags1$conf[df.tags1$cell.id %in% select.barcode] <- 1
  df.tags1$expand_conf <- rep(0, nrow(df.tags1))
  table_prop <- table(df.tags1$scRef.tag[df.tags1$conf == 1])/length(select.barcode)
  not_extend <- names(table_prop)[table_prop > 0.2]
  for (sub_cluster in unique(df.tags1$cluster.merge.id)) {
    df_sub <- df.tags1[df.tags1$cluster.merge.id == sub_cluster,]
    kind_num <- length(table(df_sub$scRef.tag))
    sub_tag <- names(sort(table(df_sub$scRef.tag), decreasing = T))[1]
    if (kind_num == 1 & sum(df_sub$conf) > 0 & !(sub_tag %in% not_extend)) {
      df.tags1[df.tags1$cluster.merge.id == sub_cluster, 'expand_conf'] <- 1
    } else {
      df.tags1[df.tags1$cluster.merge.id == sub_cluster, 'expand_conf'] <- df_sub$conf
    }
  }
  select.barcode <- df.tags1$cell.id[df.tags1$expand_conf == 1]
  select.exp <- df.exp.merge[, cell_ids %in% select.barcode]
  select.tag1 <- tag1[tag1[, 'cell_id'] %in% select.barcode, ]
  LocalRef <- generate_ref(select.exp, select.tag1,  min_cell = min_cell)
  vec.tag1 <- select.tag1[, 'tag']
  print('Cell types in local reference:')
  print(dimnames(LocalRef)[[2]])
  
  #####
  gc()
  #####
  # find local marker genes
  print('find local marker genes')
  suppressMessages(
    out.markers <-
      .find_markers_sc(
        select.exp, vec.tag1, LocalRef,
        seurat.out.group, list.localNeg,
        use_RUVseq = use_RUVseq,
        percent.high.exp = percent.high.exp,
        num_threads = num_threads
      )
  )
  local.cell.genes <- out.markers[['list.cell.genes']]
  local.cell.genes_merge <- list()
  for (one_cell in names(local.cell.genes)) {
    local.cell.genes_merge[[one_cell]] <-
      setdiff(union(list.cell.genes[[one_cell]], local.cell.genes[[one_cell]]), NA)
  }
  
  print('Second-round annotation:')
  if (method2 == 'randomforest') {
    pca_query <- as.data.frame(seurat.query@reductions$pca@cell.embeddings)
    train_pca <- pca_query[colnames(select.exp), ]
    train_pca$label <- as.factor(vec.tag1)
    library(randomForest)
    fit.forest <- randomForest(label ~ ., data = train_pca)
    pred_tags <- predict(fit.forest, pca_query)
    tag2 <- data.frame(cell_id = rownames(pca_query), tag = pred_tags)
    tag2 <- as.matrix(tag2)
  } else {
    if (!is.null(corr_use_HVGene2)) {
      num_feature <- corr_use_HVGene2/2
      train_set <- query_set[colnames(select.exp),]
      train_set$label <- vec.tag1
      HVG <- SelectGene(train_set, k = num_feature)
      similarity.in <- df.exp.merge[HVG, ]
      ref.in <- LocalRef[HVG, ]
    } else {
      similarity.in <- df.exp.merge
      ref.in <- LocalRef
    }
    if (method2 != 'multinomial') {
      out2 <- .get_cor(ref.in, similarity.in, method = method2, num_threads = num_threads)
      out2 <- t(out2)
      tag2 <- .get_tag_max(t(out2))
    } else {
      test_set <- query_set
      pred_tags <- SciBet(train_set, test_set, k=num_feature)
      tag2 <- data.frame(cell_id = rownames(query_set), tag = pred_tags)
      tag2 <- as.matrix(tag2)
    }
  }
  
  gc()
  
  list_out_2 <- .confirm_label_auc(df.exp.merge, local.cell.genes_merge, list_near_cell,
                                   tag2, num_threads = num_threads)
  df.tags2 <- list_out_2$meta.tag
  list_tags2_back <- list_out_2$list.auc.back
  cell_ids <- colnames(df.exp.merge)
  df.tags2 <- df.tags2[cell_ids, ]
  gc()
  
  
  pvalue1 <- data.frame(stringsAsFactors = F)
  for (one_cell in unique(df.tags1$scRef.tag)) {
    df_sub <- df.tags1[df.tags1$scRef.tag == one_cell, c('scRef.tag', 'AUC')]
    df_sub$tag_unassaigned <- df_sub$scRef.tag
    df_sub$tag_unassaigned[df_sub$AUC < vec.cut_1[one_cell]] <- 'Unassigned'
    pvalue1 <- rbind(pvalue1, df_sub)
  }
  names(pvalue1) <- c('scRef.tag.1', 'AUC.1', 'tag_unassaigned.1')
  df_unassigned <- pvalue1[cell_ids, ]
  
  pvalue1 <- data.frame(stringsAsFactors = F)
  for (cell_label in unique(df.tags1$scRef.tag)) {
    df_sub <- df.tags1[df.tags1$scRef.tag == cell_label,]
    df_sub$AUC_scale <- scale(df_sub$AUC)
    pvalue1 <- rbind(pvalue1, df_sub[, c('scRef.tag', 'AUC_scale')])
  }
  names(pvalue1) <- c('scRef.tag.1', 'AUC.1')
  pvalue2 <- data.frame(stringsAsFactors = F)
  for (cell_label in unique(df.tags2$scRef.tag)) {
    df_sub <- df.tags2[df.tags2$scRef.tag == cell_label,]
    df_sub$AUC_scale <- scale(df_sub$AUC)
    pvalue2 <- rbind(pvalue2, df_sub[, c('scRef.tag', 'AUC_scale')])
  }
  names(pvalue2) <- c('scRef.tag.2', 'AUC.2')
  pvalue <- merge(pvalue1, pvalue2, by = 'row.names')
  row.names(pvalue) <- pvalue$Row.names
  pvalue$Row.names <- NULL
  pvalue <- pvalue[cell_ids, ]
  mtx.tag <- as.matrix(pvalue[, c('scRef.tag.1', 'scRef.tag.2')])
  mtx.pval <- as.matrix(pvalue[, c('AUC.1', 'AUC.2')])
  mtx.rank <- apply(mtx.pval, 1, rank, ties.method = "first")
  tag.merge <-
    apply(as.array(1:dim(mtx.tag)[1]), 1, function(i) {
      mtx.tag[i, mtx.rank[2, i]]
    })
  pvalue$tag.merge <- tag.merge
  df_auc_unassigned <- merge(pvalue, df_unassigned, by = 'row.names')
  rownames(df_auc_unassigned) <- df_auc_unassigned$Row.names
  df_auc_unassigned$Row.names <- NULL
  df_auc_unassigned <- df_auc_unassigned[cell_ids, ]
  df_auc_unassigned$tag.final <- tag.merge
  if (identify_unassigned) {
    df_auc_unassigned$tag.final[df_auc_unassigned$tag_unassaigned.1 == 'Unassigned'] <- 'Unassigned'
  }
  
  if (opt_speed) {
    df.cluster <- df.dict[, c("cluster.merge.id", "cluster.level1")]
    df.cluster <- unique(df.cluster)
    df.cluster <- data.frame(cluster.id = df.cluster$cluster.level1,
                             row.names = df.cluster$cluster.merge.id,
                             stringsAsFactors = F)
  }
  
  df.tags <- merge(df_auc_unassigned, df.dict, by = 'row.names')
  row.names(df.tags) <- df.tags$Row.names
  df.tags$Row.names <- NULL
  df.tags <- df.tags[, c("tag.merge", "tag.final", "cluster.merge.id", "cluster.level1")]
  df.tags$scRef.tag <- df.tags$tag.final
  all_sub_clusters <- unique(df.tags$cluster.merge.id)
  for (sub_cluster in all_sub_clusters) {
    df_sub <- df.tags[df.tags$cluster.merge.id == sub_cluster,]
    sub_table <- sort(table(df_sub$tag.final), decreasing = T)
    if (length(sub_table) == 1) {
      next
    } else {
      if (sub_table[1]/nrow(df_sub) >= 0.6) {
        df.tags$scRef.tag[df.tags$cluster.merge.id == sub_cluster] <- names(sub_table)[1]
      }
    }
  }
  
  # if (opt_speed) {
  #     df.tags$cluster.merge.id <- row.names(df.tags)
  #     df.tags.merge <- merge(df.tags, df.dict[, c('cluster.merge.id', 'cell.id')],
  #                            by = 'cluster.merge.id')
  #     df.tags.merge$cluster.merge.id <- NULL
  #     row.names(df.tags.merge) <- df.tags.merge$cell.id
  #     df.tags.merge$cell.id <- NULL
  #     df.tags <- df.tags.merge
  # }
  cell_ids <- colnames(exp_sc_mat)
  # df.combine <- df.combine[cell_ids, ]
  df.tags <- df.tags[cell_ids, ]
  
  if (cluster_assign) {
    all_clusters <- unique(df.tags$cluster.level1)
    df.tags$cluster.tags <- df.tags$scRef.tag
    for (cluster in all_clusters) {
      df_sub <- df.tags[df.tags$cluster.level1 == cluster,]
      sub_table <- sort(table(df_sub$tag.final), decreasing = T)
      df.tags$cluster.tags[df.tags$cluster.level1 == cluster] <- names(sub_table)[1]
    }
    df.combine <- data.frame(scMAGIC.tag = df.tags$cluster.tags, row.names = rownames(df.tags))
  } else {
    df.combine <- data.frame(scMAGIC.tag = df.tags$scRef.tag, row.names = rownames(df.tags))
  }
  
  
  #####
  gc()
  #####
  time2 <- Sys.time()
  time.scMAGIC <- difftime(time2, time1, units = 'secs')
  if (simple_output) {
    output <- df.combine
  } else {
    output <- list()
    output$tag1 <- tag1
    output$out1 <- out1
    output$LocalRef <- LocalRef
    if (identify_unassigned) {
      output$pvalue1 <- df.tags1
      output$pvalue2 <- df.tags
      # output$info.cluster <- info.cluster
      if (opt_speed) {
        output$dict.cluster <- df.dict
      }
      # output$cutoff.1 <- df.cutoff.1
      # output$cutoff.neg.1 <- neg.cutoff.1
      # output$cutoff.2 <- df.cutoff.2
      output$ref.markers <- list.cell.genes
      output$local.markers <- local.cell.genes
    }
    output$final.out <- df.combine
    output$run.time <- time.scMAGIC
    if (mod == 'debug') {
      output$df.exp.merge <- df.exp.merge
      output$diff.log10Pval <- diff.log10Pval
    }
  }
  
  print('Finish!')
  
  return(output)
  
}


transformHomoloGene <- function(exp_sc_mat, inTaxID = 9606, outTaxID = 10090) {
  library(homologene)
  genes.in <- rownames(exp_sc_mat)
  res.home <- homologene(genes.in, inTax = inTaxID, outTax = outTaxID)
  res.home <- res.home[!duplicated(res.home[, 1]),]
  res.home <- res.home[!duplicated(res.home[, 2]),]
  genes.out <- res.home[, 1]
  genes.homo <- res.home[, 2]
  exp.out <- exp_sc_mat[genes.out,]
  rownames(exp.out) <- genes.homo
  return(exp.out)
}


scMAGIC_Seurat <- function(seurat.query, seurat.ref, atlas = 'MCA', corr_use_HVGene = 2000,
                           num_threads = 4) {
  exp_sc_mat <- seurat.query@assays$RNA@counts
  exp_ref_mat <- seurat.ref@assays$RNA@counts
  exp_ref_label <- seurat.ref@meta.data$celltype
  output.scMAGIC <- scMAGIC(exp_sc_mat, exp_ref_mat, exp_ref_label,
                            atlas = atlas, corr_use_HVGene1 = corr_use_HVGene,
                            corr_use_HVGene2 = corr_use_HVGene,
                            num_threads = num_threads)
  pred.scMAGIC <- output.scMAGIC$scMAGIC.tag
  seurat.query$prediction_celltype <- pred.scMAGIC
  return(seurat.query)
}

########################################
#Author: Yu Zhang
#Email: zhang_yu18@fudan.edu.cn
#######################################

.combine_tags <- function(df.tags1, df.tags2) {
  # concat reference pval and local pval
  pvalue1 <- df.tags1[, c('scRef.tag', 'pvalue')]
  names(pvalue1) <- c('scRef.tag.1', 'pvalue.1')
  pvalue2 <- df.tags2[, c('scRef.tag', 'pvalue')]
  names(pvalue2) <- c('scRef.tag.2', 'pvalue.2')
  pvalue <- merge(pvalue1, pvalue2, by = 'row.names')
  row.names(pvalue) <- pvalue$Row.names
  pvalue$Row.names <- NULL
  
  # select more confident tag
  mtx.tag <- as.matrix(pvalue[, c('scRef.tag.1', 'scRef.tag.2')])
  mtx.pval <- as.matrix(pvalue[, c('pvalue.1', 'pvalue.2')])
  mtx.rank <- apply(mtx.pval, 1, rank, ties.method = "first")
  tag.final <-
    apply(as.array(1:dim(mtx.tag)[1]), 1, function(i) {
      mtx.tag[i, mtx.rank[1, i]]
    })
  pval.final <-
    apply(as.array(1:dim(mtx.pval)[1]), 1, function(i) {
      mtx.pval[i, mtx.rank[1, i]]
    })
  tag.final <-
    data.frame(
      scRef.tag = tag.final,
      pvalue = pval.final,
      row.names = dimnames(mtx.tag)[[1]],
      stringsAsFactors = F
    )
  
  OUT <- list()
  OUT$pvalue <- pvalue
  OUT$tag.final <- tag.final
  return(OUT)
  
}


.find_markers_cell_1 <- function(cell.ref, mtx.combat, exp_ref_mat, percent.high.exp,
                                 mtx.combat.use, topN, i) {
  cell <- cell.ref[i]
  # print(cell)
  vec.cell <- mtx.combat[, paste0('Ref.', cell)]
  # vec.cell.high <- vec.cell[vec.cell > quantile(vec.cell, percent.high.exp)]
  vec.ref <- exp_ref_mat[, cell]
  vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, percent.high.exp)]
  # high expression genes
  # genes.high <- intersect(names(vec.cell.high), names(vec.ref.high))
  genes.high <- names(vec.ref.high)
  
  # function
  .getDEgeneF <- function(esetm = NULL, group = NULL, pair = FALSE,
                          block = NULL, p_adj = "fdr", fpkm = T) {
    # limma function
    if (is.null(esetm)) {
      cat(
        "esetm: gene expression matrix",
        "group: factor: \"c\"/\"d\"",
        "pair: TRUE/FALSE*",
        "block: e.g.1 2 2 1 if paired; blank if not",
        "p_adj: p.adjust, fdr* ",
        "fpkm: TRUE/FALSE*",
        sep = "\n"
      )
    } else{
      library(limma, verbose = F)
      if (pair) {
        design <- model.matrix( ~ block + group)
      } else{
        design <- model.matrix( ~ group)
      }
      suppressWarnings(fit <- lmFit(esetm, design))
      if (fpkm) {
        suppressWarnings(fit <- eBayes(fit, trend = T, robust = T))
      } else{
        suppressWarnings(fit <- eBayes(fit))
      }
      x <- topTable(fit, number = nrow(esetm), adjust.method = p_adj,
                    coef = "group2")
      x <- x[!is.na(row.names(x)), ]
      x <- x[!duplicated(row.names(x)), ]
      return(x)
    }
  }
  
  # diff in reference
  cells <- c(setdiff(cell.ref, cell), cell)
  bool.cell <- cells
  bool.cell[bool.cell != cell] <- '1'
  bool.cell[bool.cell == cell] <- '2'
  res.limma.ref <- .getDEgeneF(exp_ref_mat[genes.high, cells], bool.cell)
  res.limma.ref <- res.limma.ref[res.limma.ref$logFC > 0,]
  genes.ref <- row.names(res.limma.ref[(res.limma.ref$P.Value < 0.05),])
  if (length(genes.ref) < (3*topN)) {
    res.limma.ref <- res.limma.ref[order(res.limma.ref$P.Value),]
    genes.ref <- row.names(res.limma.ref)[1:(3*topN)]
  }
  use.genes <- intersect(genes.high, genes.ref)
  
  # diff in Atlas
  if (length(use.genes) > topN) {
    mtx.limma <- cbind(mtx.combat.use, vec.cell)
    bool.atlas.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
    res.limma.MCA <- .getDEgeneF(mtx.limma[use.genes, ], bool.atlas.cell)
    res.limma.MCA <- res.limma.MCA[res.limma.MCA$logFC > 0,]
    genes.diff <- row.names(res.limma.MCA[(res.limma.MCA$P.Value < 0.001),])
    if (length(genes.diff) < (topN)) {
      res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
      genes.diff <- row.names(res.limma.MCA)[1:(topN)]
    }
    if (length(genes.diff) > (3*topN)) {
      res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
      genes.diff <- row.names(res.limma.MCA)[1:(3*topN)]
    }
  } else {
    genes.diff <- use.genes
  }
  
  return(genes.diff)
  
}

.find_markers_1 <- function(exp_ref_mat, seurat.out.group,
                            type_ref = 'sum-counts', use_RUVseq = T,
                            base.topN = 50, percent.high.exp = 0.8, num_threads = 6) {
  library(parallel, verbose = F)
  library(Seurat, verbose = F)
  # check parameters
  if (!is.null(base.topN)) {
    base.topN <- base.topN
  } else {
    stop('Error in finding markers: provide incorrect parameters')
  }
  
  # transform count to fpm
  if (type_ref == 'sum-counts') {
    seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
    seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                scale.factor = 1e6, verbose = F)
    exp_ref_mat <- as.matrix(seurat.Ref@assays$RNA@data)
  }
  if (type_ref %in% c('fpkm', 'tpm', 'rpkm', 'bulk')) {
    exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
  }
  
  ###### regard a outgroup (e.g. MCA/HCA) as reference of DEG
  seurat.MCA <- seurat.out.group
  fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
  
  # overlap genes
  fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
  out.overlap <- get_overlap_genes(fpm.MCA, exp_ref_mat)
  fpm.MCA <- as.matrix(out.overlap$exp_sc_mat)
  exp_ref_mat <- as.matrix(out.overlap$exp_ref_mat)
  if (use_RUVseq) {
    gene_overlap <- out.overlap$gene_over
    SEG.MCA <- VariableFeatures(seurat.MCA)
    gene.constant <- intersect(gene_overlap, SEG.MCA)
  }
  # print('Number of overlapped genes:')
  # print(nrow(exp_ref_mat))
  
  cell.MCA <- dimnames(fpm.MCA)[[2]]
  cell.ref <- dimnames(exp_ref_mat)[[2]]
  # use RUVseq to remove batch effect
  # use_RUVseq=F
  mtx.in <- cbind(fpm.MCA, exp_ref_mat)
  names.mix <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
  dimnames(mtx.in)[[2]] <- names.mix
  if (use_RUVseq) {
    library(RUVSeq, verbose = F)
    seqRUVg <- RUVg(as.matrix(mtx.in), gene.constant, k=1, isLog = T)
    mtx.combat <- seqRUVg$normalizedCounts
  } else {
    mtx.combat <- mtx.in
  }
  mtx.combat <- mtx.in
  mtx.combat.use <- mtx.combat[, paste0('MCA.', cell.MCA)]
  
  
  topN <- base.topN
  cl = makeCluster(num_threads, outfile = '')
  # clusterExport(cl, '.getDEgeneF')
  RUN <- parLapply(
    cl = cl,
    1:length(cell.ref),
    .find_markers_cell_1,
    cell.ref = cell.ref, mtx.combat = mtx.combat,
    exp_ref_mat = exp_ref_mat, percent.high.exp = percent.high.exp,
    mtx.combat.use = mtx.combat.use, topN = topN
  )
  stopCluster(cl)
  names(RUN) <- cell.ref
  
  out <- list()
  out[['list.cell.genes']] <- RUN
  out[['exp_ref_mat']] <- exp_ref_mat
  return(out)
  
}


.find_markers_sc_cell_2 <- function(cell.ref, mtx.combat, LocalRef.sum, percent.high.exp,
                                    select.exp, vec.tag1, list.localNeg, mtx.combat.use, topN, i) {
  library(Seurat, verbose = F)
  cell <- cell.ref[i]
  # print(cell)
  vec.cell <- mtx.combat[, paste0('Ref.', cell)]
  # vec.cell.high <- vec.cell[vec.cell > quantile(vec.cell, percent.high.exp)]
  vec.ref <- LocalRef.sum[, cell]
  vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, percent.high.exp)]
  # high expression genes
  genes.high <- names(vec.ref.high)
  
  # function
  .getDEgeneF <- function(esetm = NULL, group = NULL, pair = FALSE,
                          block = NULL, p_adj = "fdr", fpkm = T) {
    # limma function
    if (is.null(esetm)) {
      cat(
        "esetm: gene expression matrix",
        "group: factor: \"c\"/\"d\"",
        "pair: TRUE/FALSE*",
        "block: e.g.1 2 2 1 if paired; blank if not",
        "p_adj: p.adjust, fdr* ",
        "fpkm: TRUE/FALSE*",
        sep = "\n"
      )
    } else{
      library(limma, verbose = F)
      if (pair) {
        design <- model.matrix( ~ block + group)
      } else{
        design <- model.matrix( ~ group)
      }
      suppressWarnings(fit <- lmFit(esetm, design))
      if (fpkm) {
        suppressWarnings(fit <- eBayes(fit, trend = T, robust = T))
      } else{
        suppressWarnings(fit <- eBayes(fit))
      }
      x <- topTable(fit, number = nrow(esetm), adjust.method = p_adj,
                    coef = "group2")
      x <- x[!is.na(row.names(x)), ]
      x <- x[!duplicated(row.names(x)), ]
      return(x)
    }
  }
  
  # diff in local reference and negative reference
  seurat.Ref.cell <- CreateSeuratObject(counts = select.exp, project = "Ref")
  seurat.Ref.cell <- NormalizeData(seurat.Ref.cell, verbose = F)
  seurat.Ref.cell@meta.data$original.label <- vec.tag1
  markers <- FindMarkers(seurat.Ref.cell, ident.1 = cell, group.by = 'original.label',
                         only.pos = T, features = genes.high, min.cells.group = 1,
                         min.pct = 0.3, min.diff.pct = 0.1,
                         logfc.threshold = 0.4, verbose = F)
  # markers$p_val_fdr <- p.adjust(markers$p_val, method = 'fdr')
  genes.ref <- row.names(markers[(markers$p_val_adj < 0.01),])
  if (length(genes.ref) < (3*topN)) {
    markers <- markers[order(markers$p_val),]
    genes.ref <- row.names(markers)[1:(min(3*topN, nrow(markers)))]
  }
  use.genes <- genes.ref
  # if (length(genes.ref) > (6*topN)) {
  #     markers <- markers[order(markers$p_val),]
  #     genes.ref <- row.names(res.limma.MCA)[1:(6*topN)]
  # }
  count.neg <- list.localNeg[[cell]]
  if (!is.null(count.neg)) {
    cell.exp <- select.exp[, vec.tag1 == cell]
    mat.in <- cbind(cell.exp, count.neg)
    if (sum(vec.tag1 == cell) == 1) {
      ncol.cell <- 1
    } else {
      ncol.cell <- ncol(cell.exp)
    }
    tag.in <- c(rep(cell, ncol.cell), rep('Negative', ncol(count.neg)))
    seurat.neg.cell <- CreateSeuratObject(counts = mat.in, project = "Ref")
    seurat.neg.cell <- NormalizeData(seurat.neg.cell, verbose = F)
    seurat.neg.cell@meta.data$original.label <- tag.in
    markers <- FindMarkers(seurat.neg.cell, ident.1 = cell, group.by = 'original.label',
                           only.pos = T, features = use.genes, min.cells.group = 1,
                           min.pct = 0.3, min.diff.pct = 0.1,
                           logfc.threshold = 0.4, verbose = F)
    # markers$p_val_fdr <- p.adjust(markers$p_val, method = 'fdr')
    genes.neg <- row.names(markers[(markers$p_val_adj < 0.01),])
    # if (length(genes.ref) < (2*topN)) {
    #     markers <- markers[order(markers$p_val),]
    #     genes.ref <- row.names(markers)[1:(min(2*topN, nrow(markers)))]
    # }
    if (length(genes.neg) < (2*topN)) {
      markers <- markers[order(markers$p_val),]
      genes.neg <- row.names(markers)[1:(min(2*topN, nrow(markers)))]
    }
    if (length(genes.neg) > (4*topN)) {
      markers <- markers[order(markers$p_val),]
      genes.neg <- row.names(markers)[1:(4*topN)]
    }
    use.genes <- genes.neg
  }
  
  # diff in Atlas
  if (length(use.genes) > topN) {
    mtx.limma <- cbind(mtx.combat.use, vec.cell)
    bool.atlas.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
    res.limma.MCA <- .getDEgeneF(mtx.limma[use.genes, ], bool.atlas.cell)
    res.limma.MCA <- res.limma.MCA[res.limma.MCA$logFC > 0,]
    genes.diff <- row.names(res.limma.MCA[(res.limma.MCA$P.Value < 0.001),])
    if (length(genes.diff) < (topN)) {
      res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
      genes.diff <- row.names(res.limma.MCA)[1:(topN)]
    }
    if (length(genes.diff) > (3*topN)) {
      res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
      genes.diff <- row.names(res.limma.MCA)[1:(3*topN)]
    }
  } else {
    genes.diff <- use.genes
  }
  
  return(genes.diff)
  
}


.find_markers_sc_2 <- function(select.exp, vec.tag1, LocalRef,
                               seurat.out.group, list.localNeg,
                               use_RUVseq = T, num_threads = 6,
                               base.topN = 50, percent.high.exp = 0.80) {
  library(parallel, verbose = F)
  library(Seurat, verbose = F)
  # check parameters
  if (!is.null(base.topN)) {
    base.topN <- base.topN
  } else {
    stop('Error in finding markers: provide incorrect parameters')
  }
  
  # transform count to fpm
  seurat.Ref.sum <- CreateSeuratObject(counts = LocalRef, project = "Ref")
  seurat.Ref.sum <- NormalizeData(seurat.Ref.sum, verbose = F)
  LocalRef.sum <- as.matrix(seurat.Ref.sum@assays$RNA@data)
  
  ###### regard a outgroup (e.g. MCA/HCL) as reference of DEG
  seurat.MCA <- seurat.out.group
  fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
  
  # overlap genes
  fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
  out.overlap <- get_overlap_genes(fpm.MCA, LocalRef.sum)
  fpm.MCA <- as.matrix(out.overlap$exp_sc_mat)
  LocalRef.sum <- as.matrix(out.overlap$exp_ref_mat)
  if (use_RUVseq) {
    gene_overlap <- out.overlap$gene_over
    SEG.MCA <- VariableFeatures(seurat.MCA)
    gene.constant <- intersect(gene_overlap, SEG.MCA)
  }
  # print('Number of overlapped genes:')
  # print(nrow(exp_ref_mat))
  
  cell.MCA <- dimnames(fpm.MCA)[[2]]
  cell.ref <- dimnames(LocalRef.sum)[[2]]
  # use RUVseq to remove batch effect
  mtx.in <- cbind(fpm.MCA, LocalRef.sum)
  names.mix <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
  dimnames(mtx.in)[[2]] <- names.mix
  if (use_RUVseq) {
    library(RUVSeq, verbose = F)
    seqRUVg <- RUVg(as.matrix(mtx.in), gene.constant, k=3, isLog = T)
    mtx.combat <- seqRUVg$normalizedCounts
  } else {
    mtx.combat <- mtx.in
  }
  mtx.combat.use <- mtx.combat[, paste0('MCA.', cell.MCA)]
  
  
  topN <- base.topN
  cl = makeCluster(num_threads, outfile = '')
  # clusterExport(cl, '.getDEgeneF')
  RUN <- parLapply(
    cl = cl,
    1:length(cell.ref),
    .find_markers_sc_cell_2,
    cell.ref = cell.ref, mtx.combat = mtx.combat,
    LocalRef.sum = LocalRef.sum, percent.high.exp = percent.high.exp,
    select.exp = select.exp, vec.tag1 = vec.tag1,
    list.localNeg = list.localNeg, mtx.combat.use = mtx.combat.use, topN = topN
  )
  stopCluster(cl)
  names(RUN) <- cell.ref
  
  out <- list()
  out[['list.cell.genes']] <- RUN
  out[['exp_ref_mat']] <- LocalRef.sum
  return(out)
  
}


.cluster_sc_one_out <- function(exp_sc_mat, cluster_num_pc = 75, cluster_resolution = 3) {
  # cluster
  library(Seurat, verbose = F)
  # data preparing
  seurat.unlabeled <- CreateSeuratObject(counts = exp_sc_mat)
  seurat.unlabeled <-
    NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize",
                  scale.factor = 10000, verbose = F)
  # print(seurat.unlabeled@assays$RNA@data[1:6,1:6])
  seurat.unlabeled <-
    FindVariableFeatures(seurat.unlabeled, selection.method = "vst",
                         nfeatures = 2000, verbose = F)
  seurat.unlabeled <- ScaleData(seurat.unlabeled, verbose = F)
  
  # PCA
  seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = cluster_num_pc, verbose = F)
  
  # cluster
  seurat.unlabeled <-
    FindNeighbors(seurat.unlabeled, reduction = "pca", dims = 1:cluster_num_pc,
                  nn.eps = 0.5, verbose = F)
  seurat.unlabeled <-
    FindClusters(seurat.unlabeled, resolution = cluster_resolution,
                 n.start = 10, n.iter = 100, verbose = F)
  
  out.cluster <-
    data.frame(
      cluster.id = as.character(seurat.unlabeled@meta.data$seurat_clusters),
      # original.label = seurat.unlabeled@meta.data$original.label,
      row.names = dimnames(seurat.unlabeled@assays$RNA@counts)[[2]]
    )
  return(out.cluster)
  
}


.one_confirm_label <- function(barcode, exp_sc_mat, df.tag, list.cell.genes,
                               method.test = 'wilcox.test') {
  expression.barcode <- exp_sc_mat[, barcode]
  bool.mark.gene <- rep(1, dim(exp_sc_mat)[1])
  cell <- df.tag[barcode, 'scRef.tag']
  genes.marker <- list.cell.genes[[cell]]
  bool.mark.gene[dimnames(exp_sc_mat)[[1]] %in% genes.marker] <- 2
  test.in <- cbind(expression.barcode, bool.mark.gene)
  test.in <- as.data.frame(test.in)
  names(test.in) <- c('expression.level', 'factor.mark.gene')
  test.in$factor.mark.gene <- as.factor(test.in$factor.mark.gene)
  if (method.test == 'ks.test') {
    vec.other <- test.in[test.in$factor.mark.gene == 1, 'expression.level']
    vec.marker <- test.in[test.in$factor.mark.gene == 2, 'expression.level']
    out.test <-
      ks.test(
        x = vec.other, y = vec.marker,
        data = test.in,
        alternative = 'greater'
      )
  } else {
    if (method.test == 'wilcox.test') {
      out.test <-
        wilcox.test(
          formula = expression.level ~ factor.mark.gene,
          data = test.in,
          digits.rank = 0,
          alternative = 'less'
        )
    } else {
      stop('Error: Please provide correct test method!')
    }
  }
  pvalue <- max(out.test$p.value, 1e-200)
  gc()
  return(data.frame(pvalue = pvalue, row.names = barcode))
  
}


.confirm_label <- function(exp_sc_mat, list.cell.genes, scRef.tag, num_threads = 10) {
  library(parallel, verbose = F)
  # confirm label
  exp_sc_mat <- as.matrix(exp_sc_mat)
  df.tag <- as.data.frame(scRef.tag)
  row.names(df.tag) <- df.tag$cell_id
  df.tag$cell_id <- NULL
  names(df.tag) <- 'scRef.tag'
  
  cl = makeCluster(num_threads, outfile = '')
  RUN <- parLapply(
    cl = cl,
    row.names(df.tag),
    .one_confirm_label,
    exp_sc_mat = exp_sc_mat,
    df.tag = df.tag,
    list.cell.genes = list.cell.genes
  )
  stopCluster(cl)
  df.pval = data.frame()
  for (single.pval in RUN) {
    df.pval = rbind(df.pval, single.pval)
  }
  
  meta.tag <- merge(df.tag, df.pval, by = 'row.names')
  row.names(meta.tag) <- meta.tag$Row.names
  meta.tag$Row.names <- NULL
  meta.tag$log10Pval <- -log10(meta.tag$pvalue)
  
  gc()
  
  return(meta.tag)
  
}


.cutoff_GMM_add_neg <- function(df.tags.in, num_cluster = NULL, cutoff.neg = 5,
                                cutoff.pos = 5, ceiling.cutoff = 30,
                                opt.strict = F, opt.negative = F) {
  library(mclust, verbose = F)
  cells <- as.character(unique(df.tags.in$scRef.tag))
  vec.cutoff <- c()
  vec.neg.cutoff <- c()
  for (i in 1:length(cells)) {
    cell <- cells[i]
    print(cell)
    df.sub <- df.tags.in[df.tags.in$scRef.tag == cell, ]
    if (length(unique(df.sub$log10Pval)) < 5) {
      vec.cutoff <- c(vec.cutoff, (ceiling.cutoff + cutoff.pos)/2)
      if (opt.negative) {
        neg.cutoff <- NA
        vec.neg.cutoff <- c(vec.neg.cutoff, neg.cutoff)
      }
      next()
    }
    if (is.null(num_cluster)) {
      model <- densityMclust(df.sub$log10Pval, verbose = F)
    } else {
      model <- densityMclust(df.sub$log10Pval,
                             G = min(num_cluster, round(length(unique(df.sub$log10Pval))/5)),
                             verbose = F)
    }
    cluster.mean <- model$parameters$mean[order(model$parameters$mean)]
    names(cluster.mean) <- as.character(1:length(cluster.mean))
    cluster.sd <- sqrt(model$parameters$variance$sigmasq)[order(model$parameters$mean)]
    if (length(cluster.sd) == 1) {
      cluster.sd <- rep(sqrt(model$parameters$variance$sigmasq), length(cluster.mean))
    }
    names(cluster.sd) <- names(cluster.mean)
    df.sub$classification <- model$classification
    
    # negative cutoff
    if (opt.negative) {
      cutoff.neg <- max(0, cutoff.neg-2)
      cluster.neg <- names(cluster.mean[cluster.mean < cutoff.neg])
      if (length(cluster.neg) == 0) {
        neg.cutoff <- NA
      } else {
        percent.neg <-
          sum(df.sub$log10Pval < cutoff.neg) / sum(df.sub$log10Pval < cutoff.neg + 5)
        if (percent.neg > 0.7) {
          neg.cutoff <- cutoff.neg
        } else {
          neg.cutoff <- NA
        }
      }
      vec.neg.cutoff <- c(vec.neg.cutoff, neg.cutoff)
    }
    
    
    # positive cutoff
    all.clusters <- names(cluster.mean)
    cluster.unknown <- names(cluster.mean[cluster.mean < cutoff.pos])
    if (opt.strict) {
      cluster.known <- names(cluster.mean[cluster.mean > ceiling.cutoff])
      if (sum(model$classification %in% cluster.known) < 5) {
        clusters.known <- setdiff(setdiff(all.clusters, cluster.unknown), cluster.known)
        if (length(clusters.known) > 0) {
          cluster.known <- c(clusters.known[length(clusters.known)], cluster.known)
        }
      }
    } else {
      cluster.known <- setdiff(all.clusters, cluster.unknown)
      cluster.final <- c(cluster.known[length(cluster.known)])
      if (length(cluster.known) > 1) {
        for (j in rev(1:(length(cluster.known) - 1))) {
          cluster <- cluster.known[j]
          if (cluster.mean[cluster] > ceiling.cutoff) {
            cluster.final <- c(cluster.final, cluster)
            next()
          }
          if ((cluster.mean[cluster] + 20) > cluster.mean[cluster.known[j + 1]]) {
            cluster.final <- c(cluster.final, cluster)
          } else {
            break()
          }
        }
      }
      cluster.known <- rev(cluster.final)
    }
    
    if (length(cluster.known) == 0) {
      vec.cutoff <- c(vec.cutoff, ceiling.cutoff)
      next()
    }
    if (opt.strict) {
      # sub.cutoff <- cluster.mean[cluster.known[1]]
      sub.cutoff <- median(df.sub[df.sub$classification == cluster.known[1], 'log10Pval'])
      if (sum(model$classification %in% cluster.known) < 10) {
        sub.cutoff <- min(df.sub[df.sub$classification == cluster.known[1], 'log10Pval'])
      }
    } else {
      sub.cutoff <-
        min(df.sub[df.sub$classification %in% cluster.known, 'log10Pval'])
    }
    # sub.cutoff <- min(df.sub[df.sub$classification %in% cluster.known, 'log10Pval'])
    sub.cutoff <- max(sub.cutoff, cutoff.pos)
    sub.cutoff <- min(ceiling.cutoff, sub.cutoff)
    vec.cutoff <- c(vec.cutoff, sub.cutoff)
  }
  names(vec.cutoff) <- cells
  out.cutoff <- list()
  out.cutoff$vec.cutoff <- vec.cutoff
  if (opt.negative) {
    names(vec.neg.cutoff) <- cells
    out.cutoff$vec.neg.cutoff <- vec.neg.cutoff
  }
  return(out.cutoff)
  
}


scMAGIC_atlas <- function(exp_sc_mat, exp_ref_mat, exp_ref_label = NULL,
                          type_ref = 'sum-counts', single_round = F,
                          identify_unassigned = T, atlas = 'MCA', use_RUVseq = T,
                          cluster_num_pc = 50, cluster_resolution = 3,
                          opt_speed = NULL, combine_num_cell = NULL, min_cell = NULL,
                          method1 = 'kendall', method2 = 'multinomial',
                          corr_use_HVGene1 = 2000, corr_use_HVGene2 = 2000,
                          GMM.num_component = NULL, GMM.neg_cutoff = NULL,
                          GMM.floor_cutoff = 5, GMM.ceiling_cutoff = 30,
                          threshold_recall = 0.2, num_threads = 4, simple_output = T) {
  library(parallel, verbose = F)
  # check parameters
  if (!type_ref %in% c('sc-counts', 'sum-counts', 'fpkm', 'tpm', 'rpkm')) {
    stop('Error: inexistent input of reference data format')
  }
  cutoff.1 = 'default'
  cutoff.2 = 'default'
  mod = ''
  # simple_output = T
  out.group = atlas
  if (is.null(GMM.neg_cutoff)) {
    GMM.neg_cutoff <- GMM.floor_cutoff
  }
  num.cell <- dim(exp_sc_mat)[2]
  if (is.null(opt_speed)) {
    if (num.cell > 5000) {
      opt_speed <- TRUE
    } else {
      opt_speed <- FALSE
    }
  }
  if (opt_speed) {
    if (is.null(combine_num_cell)) {
      combine_num_cell <- ceiling(num.cell/3000)
    }
    if (is.null(min_cell)) {
      min_cell <- combine_num_cell * 2
    }
  } else {
    if (is.null(min_cell)) {
      min_cell <- 1
    }
    combine_num_cell <- 1
  }
  
  time1 <- Sys.time()
  # get sum-counts format
  if (type_ref == 'sc-counts') {
    print('Sum single cell counts matrix:')
    label.in <- data.frame(cell_id = colnames(exp_ref_mat),
                           tag = as.character(exp_ref_label))
    exp_ref_mat.sum <- generate_ref(exp_ref_mat, label.in, M = 'SUM')
    exp_ref_mat <- exp_ref_mat.sum
    type_ref <- 'sum-counts'
  }
  
  # get overlap genes
  out.overlap <- get_overlap_genes(exp_sc_mat, exp_ref_mat)
  exp_sc_mat <- out.overlap$exp_sc_mat
  exp_ref_mat <- out.overlap$exp_ref_mat
  gene_over <- out.overlap$gene_over
  exp_ref_mat <- exp_ref_mat[gene_over,]
  print('Number of overlapped genes:')
  print(nrow(exp_sc_mat))
  
  seurat.out.group <- .imoprt_outgroup(out.group = out.group, use_RUVseq = use_RUVseq)
  
  if (identify_unassigned) {
    # find markers of cell types in reference
    topN = 50
    percent.high.exp = 0.8
    print('Find marker genes of cell types in reference:')
    suppressMessages(
      out.markers <-
        .find_markers_1(
          exp_ref_mat,
          seurat.out.group,
          type_ref = 'sum-counts',
          use_RUVseq = use_RUVseq,
          percent.high.exp = percent.high.exp, num_threads = num_threads
        ))
    list.cell.genes <- out.markers[['list.cell.genes']]
    genes.ref <- dimnames(out.markers[['exp_ref_mat']])[[1]]
    
    # overlap genes
    gene.overlap <- intersect(gene_over, genes.ref)
    exp_sc_mat <- exp_sc_mat[gene.overlap, ]
    exp_ref_mat <- exp_ref_mat[gene.overlap, ]
    print('Number of overlapped genes:')
    print(nrow(exp_sc_mat))
    
    # cluster analysis
    print('Start clustering :')
    df.cluster <-
      .cluster_sc_one_out(exp_sc_mat,
                          cluster_num_pc = cluster_num_pc,
                          cluster_resolution = cluster_resolution)
    print('Clustering completed!')
    
    # speed calculation
    if (opt_speed) {
      print('Speed calculation by clustering:')
      out.merge <-
        .cluster_increase_speed(exp_sc_mat, df.cluster,
                                combine_num_cell = combine_num_cell, num_threads = num_threads)
      df.dict <- out.merge$df.dict
      df.exp.merge <- out.merge$df.exp.merge
      min_cell <- ceiling(min_cell / combine_num_cell)
    } else {
      df.exp.merge <- exp_sc_mat
    }
  } else {
    df.exp.merge <- exp_sc_mat
  }
  
  # rm(exp_sc_mat)
  gc()
  df.exp.merge <- as.matrix(df.exp.merge)
  
  print('First-round annotation:')
  print(method1)
  if (!is.null(corr_use_HVGene1)) {
    HVG <- .get_high_variance_genes(exp_ref_mat, type_ref = type_ref,
                                    num.genes = corr_use_HVGene1)
    HVG = intersect(HVG, rownames(df.exp.merge))
    similarity.in <- df.exp.merge[HVG, ]
    ref.in <- exp_ref_mat[HVG, ]
  } else {
    similarity.in <- df.exp.merge
    ref.in <- exp_ref_mat
  }
  if (method1 != 'multinomial') {
    suppressMessages(
      out1 <- .get_cor(similarity.in, ref.in, method = method1, num_threads = num_threads)
    )
  } else {
    suppressMessages(
      out1 <- .get_log_p_sc_given_ref(similarity.in, ref.in, num_threads = num_threads)
    )
  }
  
  gc()
  
  tag1 <- .get_tag_max(out1)
  
  if (identify_unassigned) {
    print('Build local reference')
    suppressMessages(
      df.tags1 <- .confirm_label(df.exp.merge, list.cell.genes, tag1, num_threads = num_threads)
    )
    cell_ids <- colnames(df.exp.merge)
    df.tags1 <- df.tags1[cell_ids, ]
    
    # select cutoff.1 automatitically
    if (cutoff.1 == 'default') {
      out.cutoff <- .cutoff_GMM_add_neg(df.tags1, num_cluster = GMM.num_component,
                                        cutoff.neg = GMM.neg_cutoff,
                                        cutoff.pos = GMM.floor_cutoff,
                                        ceiling.cutoff = GMM.ceiling_cutoff,
                                        opt.strict = T, opt.negative = T)
      df.cutoff.1 <- out.cutoff$vec.cutoff
      neg.cutoff.1 <- out.cutoff$vec.neg.cutoff
      # print('Default cutoff: ')
      # print(df.cutoff.1)
      
      if (single_round) {
        df.tags <- df.tags1
        df.tags$scRef.tag.1 <- as.character(df.tags$scRef.tag)
        df.tags$scRef.tag <- df.tags$scRef.tag.1
        for (cell in names(df.cutoff.1)) {
          sub.cutoff <- df.cutoff.1[cell]
          df.tags$scRef.tag[(df.tags$scRef.tag.1 == cell) &
                              (df.tags$log10Pval < sub.cutoff)] <- 'Unassigned'
        }
        
        if (opt_speed) {
          df.cluster <- df.dict[, c("cluster.merge.id", "cluster.level1")]
          df.cluster <- unique(df.cluster)
          df.cluster <- data.frame(cluster.id = df.cluster$cluster.level1,
                                   row.names = df.cluster$cluster.merge.id,
                                   stringsAsFactors = F)
        }
        df.tags <- merge(df.tags, df.cluster, by = 'row.names')
        row.names(df.tags) <- df.tags$Row.names
        df.tags$Row.names <- NULL
        
        if (opt_speed) {
          df.tags$cluster.merge.id <- row.names(df.tags)
          df.tags.merge <- merge(df.tags, df.dict[, c('cluster.merge.id', 'cell.id')],
                                 by = 'cluster.merge.id')
          df.tags.merge$cluster.merge.id <- NULL
          row.names(df.tags.merge) <- df.tags.merge$cell.id
          df.tags.merge$cell.id <- NULL
          df.tags <- df.tags.merge
        }
        
        # recall Unassigned
        df.tags$scRef.tag.pre.recall <- df.tags$scRef.tag
        cluster.ids <- unique(df.cluster$cluster.id)
        info.cluster <- data.frame(stringsAsFactors = F)
        for (cluster.id in cluster.ids) {
          sub.tag.cluster <- df.tags[df.tags$cluster.id == cluster.id,]
          table.tag <- table(sub.tag.cluster$scRef.tag.1)
          table.tag <- table.tag[order(table.tag, decreasing = T)]
          main.cell <- names(table.tag[1])
          percent.main.cell <- table.tag[1] / dim(sub.tag.cluster)[1]
          num.Unassigned <- sum(sub.tag.cluster$scRef.tag.pre.recall == 'Unassigned')
          percent.Unassigned <- num.Unassigned / nrow(sub.tag.cluster)
          if (percent.Unassigned < threshold_recall) {
            if (percent.main.cell > (1- threshold_recall - 0.05)) {
              df.tags[(df.tags$cluster.id == cluster.id) &
                        (df.tags$scRef.tag == 'Unassigned'), 'scRef.tag'] <-
                rep(main.cell, num.Unassigned)
            } else {
              df.tags[df.tags$cluster.id == cluster.id, 'scRef.tag'] <-
                df.tags[df.tags$cluster.id == cluster.id, 'scRef.tag.1']
            }
          }
          info.cluster <- rbind(
            info.cluster,
            data.frame(
              cluster.id = cluster.id,
              percent.Unassigned = percent.Unassigned,
              main.cell = main.cell,
              percent.main.cell = percent.main.cell,
              stringsAsFactors = F
            )
          )
        }
        
        df.combine <- df.tags[, c("scRef.tag", "log10Pval")]
        cell_ids <- colnames(exp_sc_mat)
        df.combine <- df.combine[cell_ids, ]
        
        gc()
        
        time2 <- Sys.time()
        time.scRef <- difftime(time2, time1, units = 'secs')
        output <- list()
        output$tag1 <- tag1
        output$out1 <- out1
        output$combine.out <- df.tags
        output$info.cluster <- info.cluster
        if (opt_speed) {
          output$dict.cluster <- df.dict
        }
        output$ref.markers <- list.cell.genes
        output$final.out <- df.combine
        output$run.time <- time.scRef
        
        print('Finish!')
        
        return(output)
        
      }
      
      select.barcode <- c()
      for (cell in names(df.cutoff.1)) {
        sub.cutoff <- df.cutoff.1[cell]
        sub.select <- df.tags1[df.tags1$scRef.tag == cell, ]
        sub.select <- sub.select[sub.select$log10Pval >= sub.cutoff, ]
        select.barcode <- c(select.barcode, row.names(sub.select))
      }
      list.localNeg <- list()
      for (cell in names(neg.cutoff.1)) {
        # print(cell)
        sub.cutoff <- neg.cutoff.1[cell]
        if (is.na(sub.cutoff)) {
          next()
        }
        sub.select <- df.tags1[df.tags1$scRef.tag == cell, ]
        sub.select <- sub.select[sub.select$log10Pval <= sub.cutoff, ]
        neg.barcode <- row.names(sub.select)
        if (length(neg.barcode) < 3) {
          next()
        }
        neg.exp <- df.exp.merge[, cell_ids %in% neg.barcode]
        list.localNeg[[cell]] <- neg.exp
      }
    } else {
      # cutoff.1 <- 1e-20
      # select.df.tags <- df.tags1[df.tags1$log10Pval >= cutoff.1, ]
      # select.barcode <- row.names(select.df.tags)
    }
    select.exp <- df.exp.merge[, cell_ids %in% select.barcode]
    select.tag1 <- tag1[tag1[, 'cell_id'] %in% select.barcode, ]
    LocalRef <- generate_ref(select.exp, select.tag1,  min_cell = min_cell)
    vec.tag1 <- select.tag1[, 'tag']
    
  } else {
    if (single_round) {
      df.combine <- as.data.frame(tag1)
      row.names(df.combine) <- df.combine$cell_id
      df.combine$cell_id <- NULL
      names(df.combine) <- 'scRef.tag'
      cell_ids <- colnames(exp_sc_mat)
      df.combine <- data.frame(scRef.tag = df.combine[cell_ids, ], row.names = cell_ids)
      
      time2 <- Sys.time()
      time.scRef <- difftime(time2, time1, units = 'secs')
      output <- list()
      output$tag1 <- tag1
      output$out1 <- out1
      output$final.out <- df.combine
      output$run.time <- time.scRef
      
      print('Finish!')
      
      return(output)
      
    }
    LocalRef <- generate_ref(df.exp.merge, tag1, min_cell = min_cell)
  }
  print('Cell types in local reference:')
  print(dimnames(LocalRef)[[2]])
  
  #####
  gc()
  #####
  if (identify_unassigned) {
    # find local marker genes
    print('find local marker genes')
    suppressMessages(
      out.markers <-
        .find_markers_sc_2(
          select.exp, vec.tag1, LocalRef,
          seurat.out.group, list.localNeg,
          use_RUVseq = use_RUVseq,
          percent.high.exp = percent.high.exp,
          num_threads = num_threads
        )
    )
    local.cell.genes <- out.markers[['list.cell.genes']]
  }
  
  print('Second-round annotation:')
  print(method2)
  if (!is.null(corr_use_HVGene2)) {
    HVG <- .get_high_variance_genes(LocalRef, num.genes = corr_use_HVGene2)
    similarity.in <- df.exp.merge[HVG, ]
    ref.in <- LocalRef[HVG, ]
  } else {
    similarity.in <- df.exp.merge
    ref.in <- LocalRef
  }
  if (method2 != 'multinomial') {
    suppressMessages(
      out2 <- .get_cor(similarity.in, ref.in, method = method2, num_threads = num_threads)
    )
  } else {
    suppressMessages(
      out2 <- .get_log_p_sc_given_ref(similarity.in, ref.in, num_threads = num_threads)
    )
  }
  tag2 <- .get_tag_max(out2)
  gc()
  
  if (identify_unassigned) {
    suppressMessages(
      df.tags2 <- .confirm_label(df.exp.merge, local.cell.genes, tag2, num_threads = num_threads)
    )
    cell_ids <- colnames(df.exp.merge)
    df.tags2 <- df.tags2[cell_ids, ]
    gc()
    
    # combine reference pval and local pval
    del.cells <- setdiff(colnames(exp_ref_mat), colnames(LocalRef))
    df.tags1$pvalue[df.tags1$scRef.tag %in% del.cells] <- 1
    out.combine <- .combine_tags(df.tags1, df.tags2)
    tag.final <- out.combine$tag.final
    # tag.final <- df.tags2
    pvalue <- out.combine$pvalue
    
    # modify tags and combine pval
    df.tags <- merge(pvalue, tag.final, by = 'row.names')
    row.names(df.tags) <- df.tags$Row.names
    df.tags$Row.names <- NULL
    
    df.tags$log10Pval <- -log10(df.tags$pvalue)
    
    if (opt_speed) {
      df.cluster <- df.dict[, c("cluster.merge.id", "cluster.level1")]
      df.cluster <- unique(df.cluster)
      df.cluster <- data.frame(cluster.id = df.cluster$cluster.level1,
                               row.names = df.cluster$cluster.merge.id,
                               stringsAsFactors = F)
    }
    
    # select cutoff.2 automatitically
    if (cutoff.2 == 'default') {
      diff.log10Pval <- df.tags2[select.barcode, 'log10Pval'] -
        df.tags1[select.barcode, 'log10Pval']
      min.diff <- max(0, boxplot.stats(diff.log10Pval, coef = 1.5)$stats[1])
      out.cutoff.2 <- .cutoff_GMM_add_neg(df.tags, num_cluster = GMM.num_component,
                                          cutoff.neg = GMM.neg_cutoff,
                                          cutoff.pos = GMM.floor_cutoff + min.diff,
                                          ceiling.cutoff = GMM.ceiling_cutoff + min.diff)
      df.cutoff.2 <- out.cutoff.2$vec.cutoff
      # print('Default cutoff: ')
      # print(df.cutoff.2)
      df.tags$scRef.tag.12 <- as.character(df.tags$scRef.tag)
      df.tags$scRef.tag <- df.tags$scRef.tag.12
      for (cell in names(df.cutoff.2)) {
        sub.cutoff <- df.cutoff.2[cell]
        df.tags$scRef.tag[(df.tags$scRef.tag.12 == cell) &
                            (df.tags$log10Pval < sub.cutoff)] <- 'Unassigned'
      }
    } else {
      # print('Cutoff: ')
      # print(cutoff.2)
      # df.tags$scRef.tag.12 <- as.character(df.tags$scRef.tag.12)
      # df.tags$scRef.tag <- df.tags$scRef.tag.12
      # df.tags$scRef.tag[df.tags$log10Pval < cutoff.2] <- 'Unassigned'
    }
    
    df.tags <- merge(df.tags, df.cluster, by = 'row.names')
    row.names(df.tags) <- df.tags$Row.names
    df.tags$Row.names <- NULL
    
    if (opt_speed) {
      df.tags$cluster.merge.id <- row.names(df.tags)
      df.tags.merge <- merge(df.tags, df.dict[, c('cluster.merge.id', 'cell.id')],
                             by = 'cluster.merge.id')
      df.tags.merge$cluster.merge.id <- NULL
      row.names(df.tags.merge) <- df.tags.merge$cell.id
      df.tags.merge$cell.id <- NULL
      df.tags <- df.tags.merge
    }
    
    # recall Unassigned and other cell types
    df.tags$scRef.tag.pre.recall <- df.tags$scRef.tag
    cluster.ids <- unique(df.cluster$cluster.id)
    info.cluster <- data.frame(stringsAsFactors = F)
    for (cluster.id in cluster.ids) {
      sub.tag.cluster <- df.tags[df.tags$cluster.id == cluster.id,]
      num.cell <- nrow(sub.tag.cluster)
      table.tag <- table(sub.tag.cluster$scRef.tag.pre.recall)
      table.tag <- table.tag[order(table.tag, decreasing = T)]
      main.cell <- names(table.tag[1])
      percent.main.cell <- table.tag[1] / dim(sub.tag.cluster)[1]
      if ((main.cell == 'Unassigned') & (length(table.tag) > 1)) {
        main.cell <- names(table.tag[2])
        percent.main.cell <- table.tag[2] / dim(sub.tag.cluster)[1]
      }
      num.Unassigned <- sum(sub.tag.cluster$scRef.tag.pre.recall == 'Unassigned')
      percent.Unassigned <- num.Unassigned / num.cell
      num.other <- sum(sub.tag.cluster$scRef.tag.pre.recall != main.cell)
      percent.other <- num.other / num.cell
      cell.other <- paste(setdiff(names(table.tag), c(main.cell, 'Unassigned')),
                          collapse = '|')
      if (percent.Unassigned < threshold_recall) {
        if (percent.main.cell > (1 - threshold_recall - 0.2)) {
          df.tags[(df.tags$cluster.id == cluster.id) &
                    (df.tags$scRef.tag.pre.recall == 'Unassigned'), 'scRef.tag'] <-
            rep(main.cell, num.Unassigned)
        } else {
          df.tags[df.tags$cluster.id == cluster.id, 'scRef.tag'] <-
            df.tags[df.tags$cluster.id == cluster.id, 'scRef.tag.pre.recall']
        }
      }
      info.cluster <- rbind(
        info.cluster,
        data.frame(
          cluster.id = cluster.id,
          num.cell = num.cell,
          main.cell = main.cell,
          percent.main.cell = percent.main.cell,
          percent.other = percent.other,
          percent.Unassigned = percent.Unassigned,
          cell.other = cell.other,
          stringsAsFactors = F
        )
      )
    }
    
    df.combine <- df.tags[, c("scRef.tag", "log10Pval")]
    names(df.combine) <- c("scMAGIC.tag", "ConfidenceScore")
    cell_ids <- colnames(exp_sc_mat)
    df.combine <- df.combine[cell_ids, ]
    df.tags <- df.tags[cell_ids, ]
    
  } else {
    df.combine <- as.data.frame(tag2)
    row.names(df.combine) <- df.combine$cell_id
    df.combine$cell_id <- NULL
    names(df.combine) <- 'scMAGIC.tag'
    cell_ids <- colnames(exp_sc_mat)
    df.combine <- data.frame(scMAGIC.tag = df.combine[cell_ids, ], row.names = cell_ids)
  }
  
  #####
  gc()
  #####
  time2 <- Sys.time()
  time.scMAGIC <- difftime(time2, time1, units = 'secs')
  if (simple_output) {
    output <- df.combine
  } else {
    output <- list()
    output$tag1 <- tag1
    output$out1 <- out1
    output$LocalRef <- LocalRef
    if (identify_unassigned) {
      output$pvalue1 <- df.tags1
      output$pvalue2 <- df.tags
      output$info.cluster <- info.cluster
      if (opt_speed) {
        output$dict.cluster <- df.dict
      }
      output$cutoff.1 <- df.cutoff.1
      output$cutoff.neg.1 <- neg.cutoff.1
      output$cutoff.2 <- df.cutoff.2
      output$ref.markers <- list.cell.genes
      output$local.markers <- local.cell.genes
    }
    output$final.out <- df.combine
    output$run.time <- time.scMAGIC
    if (mod == 'debug') {
      output$df.exp.merge <- df.exp.merge
      output$diff.log10Pval <- diff.log10Pval
    }
  }
  
  print('Finish!')
  
  return(output)
  
}


transformHomoloGene <- function(exp_sc_mat, inTaxID = 9606, outTaxID = 10090) {
  library(homologene)
  genes.in <- rownames(exp_sc_mat)
  res.home <- homologene(genes.in, inTax = inTaxID, outTax = outTaxID)
  res.home <- res.home[!duplicated(res.home[, 1]),]
  res.home <- res.home[!duplicated(res.home[, 2]),]
  genes.out <- res.home[, 1]
  genes.homo <- res.home[, 2]
  exp.out <- exp_sc_mat[genes.out,]
  rownames(exp.out) <- genes.homo
  return(exp.out)
}

# #single cell data
# cortex_sc = Read10X("scRNAseq/filtered_feature_bc_matrix/")
# cortex_sc = CreateSeuratObject(counts = cortex_sc)
# cortex_sc
# 
# #spatial data
# anterior = Load10X_Spatial("spatialData/", filename = "filtered_feature_bc_matrix.h5")
# anterior

#extract count data from seurat object
# counts = as.matrix(cortex_sc@assays$RNA@counts)
# 
# data = data("HCL_ref")

#function for celltype identification
scm = function(single_cell_data, ref_data, c){
  pdf("scMAGIC.pdf", width = 100, height = 100)
  output.scMAGIC = scMAGIC_atlas(single_cell_data, ref_data,
                                 atlas = 'HCL',
                                 type_ref = 'sum-counts', use_RUVseq = F,
                                 min_cell = 5, num_threads = c)
  conf = output.scMAGIC$ConfidenceScore
  celltype = output.scMAGIC$scMAGIC.tag
  dfMAGIC = cbind.data.frame(colnames(single_cell_data), celltype)
  dfMAGIC = cbind.data.frame(dfMAGIC, conf)
  rownames(dfMAGIC) = dfMAGIC[,1]
  colnames(dfMAGIC) = c("ID", "scMAGIC", "conf")
  dev.off()
  
  return(dfMAGIC)
}


