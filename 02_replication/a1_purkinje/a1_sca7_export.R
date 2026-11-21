# A1 SCA7 step; run on the cluster (needs >= 16 GB RAM, Seurat + Matrix).
# GSE269430 (SCA7-266Q, Purkinje-enriched, MULTI-seq multiplexed): the two Seurat objects, 5 weeks and 8.5 weeks.
#   inspect: write META_SUMMARY.txt for each object (assays, layers, every meta.data column with its class, number of distinct
#            values and, when fewer than 60, the counts per value); nothing else.
#   export:  raw counts of the nuclei whose <cluster_col> is a Purkinje cluster in that object's label file, summed per
#            <sample_col>; genotype from <genotype_col> (a value matching SCA7 or 266Q = KI, else WT). Writes COUNTS.tsv
#            (genes x samples), SAMPLES.tsv and SHA256SUMS into a NEW <out>/ (the limma step reads it as pseudobulk_sca7/).
# v02: the objects carry no cell-type column, so the Purkinje nuclei come from the
#   CLUSTER_LABELS_<age>.tsv (columns cluster, is_purkinje; one file per object, under the pre-specified marker rule)
#   instead of <celltype_col> == <purkinje_label>, and SAMPLES.tsv records each label file's sha256; no other change.
# v03: the two objects have different gene sets, so COUNTS.tsv keeps the intersection
#   of gene names (in the 5 wk object's order); the genes dropped from each object go to DROPPED_GENES.tsv; the nine genes
#   must all be in the intersection; no other change.
# Usage (one object at a time is loaded and freed):
#   Rscript a1_sca7_export.R inspect <out_dir> <rds_5wk> <rds_8wk>
#   Rscript a1_sca7_export.R export  <out_dir> <rds_5wk> <rds_8wk> <sample_col> <genotype_col> <cluster_col> <labels_5wk> <labels_8wk>
suppressPackageStartupMessages({library(Seurat); library(Matrix)})
a <- commandArgs(trailingOnly = TRUE)
mode <- a[1]; out <- a[2]; rds <- c('5' = a[3], '8.5' = a[4])
if (dir.exists(out)) stop('output folder already exists')
sha <- function(p) strsplit(system2('sha256sum', shQuote(p), stdout = TRUE), ' ')[[1]][1]
counts_of <- function(obj) {
  assay <- DefaultAssay(obj)
  if (packageVersion('SeuratObject') >= '5.0.0') LayerData(obj, assay = assay, layer = 'counts') else GetAssayData(obj, assay = assay, slot = 'counts')
}
dir.create(out)
if (mode == 'inspect') {
  con <- file(file.path(out, 'META_SUMMARY.txt'), open = 'wx')
  for (age in names(rds)) {
    obj <- readRDS(rds[[age]]); md <- obj@meta.data
    writeLines(c(paste0('== ', basename(rds[[age]]), ' (', age, ' weeks): ', ncol(obj), ' nuclei x ', nrow(obj), ' features'),
                 paste('assays:', paste(Assays(obj), collapse = ', '), '| default:', DefaultAssay(obj),
                       '| SeuratObject', as.character(packageVersion('SeuratObject')))), con)
    for (cn in names(md)) {
      v <- md[[cn]]; u <- length(unique(v))
      writeLines(paste0('  ', cn, ' [', class(v)[1], ', ', u, ' distinct]'), con)
      if (u < 60) { tb <- sort(table(v, useNA = 'ifany'), decreasing = TRUE); writeLines(paste0('      ', names(tb), ': ', tb), con) }
    }
    rm(obj, md); invisible(gc())
  }
  close(con)
} else if (mode == 'export') {
  scol <- a[5]; gcol <- a[6]; ccol <- a[7]; lab <- c('5' = a[8], '8.5' = a[9])
  pb <- list(); meta <- list(); gn <- list()
  for (age in names(rds)) {
    obj <- readRDS(rds[[age]]); md <- obj@meta.data
    for (cn in c(scol, gcol, ccol)) if (!cn %in% names(md)) stop('missing meta.data column: ', cn)
    cnt <- counts_of(obj); gn[[age]] <- rownames(cnt)
    if (any(cnt@x != round(cnt@x))) stop('counts layer is not integer: ', basename(rds[[age]]))
    cl <- read.delim(lab[[age]], colClasses = 'character')
    if (!all(c('cluster', 'is_purkinje') %in% names(cl)) || anyDuplicated(cl$cluster)) stop('bad label file: ', lab[[age]])
    if (!all(as.character(unique(md[[ccol]])) %in% cl$cluster)) stop('a ', ccol, ' value has no label in ', lab[[age]])
    yes <- toupper(cl$is_purkinje) %in% c('TRUE', 'T', 'YES', '1')
    if (!all(yes | toupper(cl$is_purkinje) %in% c('FALSE', 'F', 'NO', '0'))) stop('is_purkinje not yes/no: ', lab[[age]])
    pc <- as.character(md[[ccol]]) %in% cl$cluster[yes]
    s <- as.character(md[[scol]])
    for (smp in sort(unique(s[pc]))) {
      j <- which(pc & s == smp); g <- unique(as.character(md[[gcol]][j]))
      if (length(g) != 1) stop('sample with more than one genotype: ', smp)
      id <- paste0(age, 'w_', smp)
      pb[[id]] <- Matrix::rowSums(cnt[, j, drop = FALSE])
      meta[[id]] <- data.frame(sample = id, genotype = if (grepl('SCA7|266Q', g, ignore.case = TRUE)) 'KI' else 'WT', genotype_raw = g,
                               age_weeks = as.numeric(age), n_nuclei = sum(s == smp), n_purkinje = length(j), kept = TRUE,
                               labels_sha256 = sha(lab[[age]]))
    }
    rm(obj, md, cnt); invisible(gc())
  }
  stopifnot(!anyDuplicated(gn[['5']]), !anyDuplicated(gn[['8.5']]))
  genes <- intersect(gn[['5']], gn[['8.5']])
  stopifnot(c('Kcnma1', 'Kcnc3', 'Kcnd3', 'Cacna1a', 'Cacna1g', 'Itpr1', 'Atp2a3', 'Kcnmb4', 'Cacna2d2') %in% genes)
  dropped <- do.call(rbind, lapply(names(gn), function(age) { d <- setdiff(gn[[age]], genes)
    data.frame(age_weeks = rep(as.numeric(age), length(d)), gene = d) }))
  con <- file(file.path(out, 'DROPPED_GENES.tsv'), open = 'wx'); write.table(dropped, con, sep = '\t', quote = FALSE, row.names = FALSE); close(con)
  m <- do.call(cbind, lapply(pb, function(x) x[genes])); rownames(m) <- genes
  con <- file(file.path(out, 'COUNTS.tsv'), open = 'wx'); write.table(m, con, sep = '\t', quote = FALSE, col.names = NA); close(con)
  sm <- do.call(rbind, meta); rownames(sm) <- sm$sample
  con <- file(file.path(out, 'SAMPLES.tsv'), open = 'wx'); write.table(sm[, -1], con, sep = '\t', quote = FALSE, col.names = NA); close(con)
  print(sm); cat("genes kept", length(genes), "| dropped 5 wk", sum(dropped$age_weeks == 5), "| dropped 8.5 wk", sum(dropped$age_weeks == 8.5), "\n")
} else stop('mode must be inspect or export')
con <- file(file.path(out, 'SHA256SUMS'), open = 'wx')
for (f in sort(list.files(out))) if (f != 'SHA256SUMS') writeLines(paste0(sha(file.path(out, f)), '  ', f), con)
close(con)
writeLines(c(paste('script', sha(normalizePath(sub('--file=', '', grep('--file=', commandArgs(FALSE), value = TRUE))))), readLines(file.path(out, 'SHA256SUMS'))))

#%% debug
