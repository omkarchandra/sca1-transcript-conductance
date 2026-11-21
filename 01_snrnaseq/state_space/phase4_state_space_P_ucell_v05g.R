#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!(length(args) %in% c(8L, 9L, 10L))) {
  stop(paste("usage:", commandArgs()[1], "counts.mtx genes.txt features.tsv chunk.nuclei.tsv random_dir|- output_prefix maxRanks chunk_size"))
}
counts_path <- args[[1]]
genes_path <- args[[2]]
features_path <- args[[3]]
metadata_path <- args[[4]]
random_dir <- args[[5]]
out_prefix <- args[[6]]
max_ranks <- as.integer(strsplit(args[[7]], ",", fixed = TRUE)[[1]])
chunk_size <- as.integer(args[[8]])
mode <- if (length(args) >= 9L) args[[9]] else "standard"
base_scores_path <- if (length(args) >= 10L) args[[10]] else "-"

suppressPackageStartupMessages({
  library(Matrix)
  library(UCell)
})
stopifnot(as.character(packageVersion("UCell")) == "2.6.2")
stopifnot(chunk_size == 256L, length(max_ranks) >= 1L, all(max_ranks %in% c(1000L, 1500L, 3000L)))

read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  fields <- strsplit(lines, "\t", fixed = TRUE)
  ids <- vapply(fields, `[[`, character(1), 1L)
  if (anyDuplicated(ids)) stop("duplicate GMT IDs in ", path)
  signatures <- lapply(fields, function(x) {
    if (length(x) < 3L) stop("invalid GMT row in ", path)
    unique(x[-c(1L, 2L)])
  })
  names(signatures) <- ids
  signatures
}

genes <- readLines(genes_path, warn = FALSE, encoding = "UTF-8")
if (!length(genes) || anyDuplicated(genes) || any(!nzchar(genes))) stop("invalid gene-symbol axis")
metadata <- read.delim(metadata_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("nucleus_id", "genotype", "age_weeks") %in% names(metadata))) stop("chunk metadata columns absent")
if (anyDuplicated(metadata$nucleus_id)) stop("duplicate nucleus IDs in chunk")

# MatrixMarket is written cells x genes by the Python H5AD streamer. UCell receives
# its native genes x cells orientation, using raw nonnegative integer counts.
counts <- Matrix::t(Matrix::readMM(counts_path))
counts <- methods::as(counts, "dgCMatrix")
if (!identical(as.integer(dim(counts)), c(length(genes), nrow(metadata)))) stop("count dimensions do not match gene/nucleus axes")
rownames(counts) <- genes
colnames(counts) <- metadata$nucleus_id

feature_table <- read.delim(features_path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("feature_id", "matched_gene_symbols") %in% names(feature_table))) stop("matched feature columns absent")
features <- strsplit(feature_table$matched_gene_symbols, ";", fixed = TRUE)
names(features) <- feature_table$feature_id
loo_mode <- identical(mode, "loo")
mouse_summary_mode <- identical(mode, "mouse_summary")
if (!mode %in% c("standard", "loo", "mouse_summary")) stop("unknown scoring mode: ", mode)
if (loo_mode && !all(grepl("^LOO__", names(features)))) stop("LOO mode requires LOO feature IDs")
if (anyDuplicated(names(features)) || (loo_mode && any(lengths(features) < 1L | lengths(features) > 199L)) ||
    (!loo_mode && any(lengths(features) < 15L | lengths(features) > 200L))) stop("invalid matched gene-set sizes")
if (any(!vapply(features, function(x) all(x %in% rownames(counts)), logical(1)))) stop("a matched set contains a symbol outside the matrix")

score_at_rank <- function(ranks, sigs) {
  value <- UCell::ScoreSignatures_UCell(precalc.ranks = ranks, features = sigs,
                                        chunk.size = chunk_size, ncores = 1, name = "")
  if (!is.matrix(value) || nrow(value) != ncol(counts)) stop("UCell returned an unexpected score matrix")
  if (any(!is.finite(value))) stop("UCell returned non-finite scores")
  value
}

ranks_by_max <- lapply(max_ranks, function(k) {
  UCell::StoreRankings_UCell(counts, maxRank = k, chunk.size = chunk_size,
                             ncores = 1, ties.method = "average", force.gc = FALSE)
})
names(ranks_by_max) <- as.character(max_ranks)
if (loo_mode) {
  if (!identical(max_ranks, 1500L) || identical(base_scores_path, "-")) stop("LOO mode requires maxRank 1500 and base scores")
  required <- c("base_feature_id", "base_feature_index", "omitted_gene")
  if (!all(required %in% names(feature_table))) stop("LOO feature table lacks base mapping columns")
  base_ids <- unique(feature_table$base_feature_id)
  n_base <- length(base_ids)
  base_raw <- readBin(base_scores_path, what = "double", n = ncol(counts) * n_base, size = 8L, endian = "little")
  if (length(base_raw) != ncol(counts) * n_base) stop("base score block has the wrong byte length")
  base_matrix <- matrix(base_raw, nrow = ncol(counts), ncol = n_base)
  wt_rows <- metadata$genotype == "WT" & as.integer(metadata$age_weeks) %in% c(5L, 18L, 24L, 30L)
  if (!any(wt_rows)) stop("LOO effects require WT discovery nuclei in this chunk")
  effects <- numeric(length(features))
  for (base_id in base_ids) {
    idx <- which(feature_table$base_feature_id == base_id)
    base_index <- unique(as.integer(feature_table$base_feature_index[idx]))
    if (length(base_index) != 1L || base_index < 0L || base_index >= n_base) stop("invalid base feature index")
    values <- score_at_rank(ranks_by_max[["1500"]], features[idx])
    base_score <- base_matrix[, base_index + 1L]
    effects[idx] <- colSums(abs(sweep(values, 1L, base_score, FUN = "-"))[wt_rows, , drop = FALSE])
  }
  con <- file(paste0(out_prefix, ".loo_abs_delta_sum.bin"), open = "wb")
  writeBin(as.double(effects), con, size = 8L, endian = "little")
  close(con)
} else if (mouse_summary_mode) {
  if (!identical(max_ranks, 1500L) || !identical(random_dir, "-")) {
    stop("mouse_summary mode requires maxRank 1500 and no random-control directory")
  }
  if (!all(c("mouse_id", "age_weeks") %in% names(metadata))) {
    stop("mouse_summary mode requires mouse_id and age_weeks metadata")
  }
  observed <- score_at_rank(ranks_by_max[["1500"]], features)
  mouse_ids <- as.character(metadata$mouse_id)
  group_ids <- unique(mouse_ids)
  group_rows <- lapply(group_ids, function(mouse_id) which(mouse_ids == mouse_id))
  sums <- do.call(rbind, lapply(group_rows, function(rows) colSums(observed[rows, , drop = FALSE])))
  ages <- vapply(group_rows, function(rows) {
    values <- unique(as.integer(metadata$age_weeks[rows]))
    if (length(values) != 1L) stop("mouse_summary chunk has inconsistent age for a mouse")
    values[[1]]
  }, integer(1))
  summary_metadata <- data.frame(mouse_id=group_ids, age_weeks=ages,
                                 n_nuclei=vapply(group_rows, length, integer(1)),
                                 stringsAsFactors=FALSE)
  write.table(summary_metadata, paste0(out_prefix, ".mouse_summary_metadata.tsv"),
              sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  con <- file(paste0(out_prefix, ".mouse_summary.bin"), open="wb")
  writeBin(as.double(sums), con, size=8L, endian="little")
  close(con)
} else {
  observed <- lapply(max_ranks, function(k) score_at_rank(ranks_by_max[[as.character(k)]], features))
  names(observed) <- as.character(max_ranks)
  for (k in max_ranks) {
    # Only primary observed scores leave R; sensitivity ranks are retained as
    # random-share counts and do not need duplicate score binaries.
    if (k == 1500L) {
      path <- paste0(out_prefix, ".score_", k, ".bin")
      con <- file(path, open = "wb")
      writeBin(as.double(observed[[as.character(k)]]), con, size = 8L, endian = "little")
      close(con)
    }
  }

  if (!identical(random_dir, "-")) {
    files <- if (dir.exists(random_dir)) {
      sort(list.files(random_dir, pattern = "^batch_[0-9]+[.]gmt$", full.names = TRUE))
    } else if (file.exists(random_dir)) {
      random_dir
    } else {
      character()
    }
    if (!length(files)) stop("random-control GMT batch directory/file is empty or absent")
    wt_discovery <- metadata$genotype == "WT" & as.integer(metadata$age_weeks) %in% c(5L, 18L, 24L, 30L)
    if (!any(wt_discovery)) stop("random controls requested for a chunk with no WT discovery nuclei")
    for (k in max_ranks) {
      rank_name <- as.character(k)
      share_counts <- matrix(0L, nrow = ncol(counts), ncol = length(features),
                             dimnames = list(colnames(counts), names(features)))
      for (batch_path in files) {
        random_features <- read_gmt(batch_path)
        parents <- sub("_R[0-9]{3}$", "", names(random_features))
        if (any(parents == names(random_features))) stop("random feature ID lacks a parent-set suffix")
        scores <- score_at_rank(ranks_by_max[[rank_name]], random_features)
        parent_order <- unique(parents)
        for (j in seq_along(parent_order)) {
          parent <- parent_order[[j]]
          random_cols <- which(parents == parent)
          if (length(random_cols) != 200L) stop("random-set count is not exactly 200 for ", parent)
          if (!parent %in% colnames(observed[[rank_name]])) stop("random parent is absent from observed UCell matrix")
          feature_column <- match(parent, names(features))
          if (is.na(feature_column)) stop("random parent is absent from the feature table: ", parent)
          # §27.1: per nucleus, retain the full share of 200 controls beaten by the observed score.
          null_scores <- scores[, random_cols, drop = FALSE]
          observed_score <- observed[[rank_name]][, parent]
          share_counts[, feature_column] <- rowSums(observed_score > null_scores)
        }
      }
      con <- file(paste0(out_prefix, ".sharecount_", k, ".bin"), open = "wb")
      writeBin(as.raw(as.integer(share_counts)), con, size = 1L)
      close(con)
    }
  }
}
cat("UCELL_VERSION=", as.character(packageVersion("UCell")),
    "\nN_CELLS=", ncol(counts), "\nN_FEATURES=", length(features),
    "\nLIBPATHS=", paste(.libPaths(), collapse = ";"), "\n", sep = "")
