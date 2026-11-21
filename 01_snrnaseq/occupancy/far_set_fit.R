# Far-set KI/WT fits, Purkinje cells .  EXPLORATORY; KI nuclei selected on the genotype lanes; ratios biased upward by construction.
# Called by far_set.py.  Same structure as the S08 fallback per-age fit (${REPO_ROOT}/01_snrnaseq/s08_limma/s08_limma_real.R): TMM + expression filter ONCE on the whole 24-mouse population matrix,
# then voom(frozen effective library sizes) -> lmFit -> contrasts.fit -> eBayes(robust=TRUE) -> topTable(confint=0.95), pooled literal 8-column ~0+genotype_age design and per-age 3v3 fits (keep.lib.sizes).
# make_design / contrast_vec / fit_tables below are lines 102-122 of the S08 script, byte for byte.
suppressPackageStartupMessages({library(limma); library(edgeR)})
args <- commandArgs(trailingOnly = TRUE); opt <- list()
for (a in args) { kv <- strsplit(sub('^--', '', a), '=', fixed = TRUE)[[1]]; opt[[kv[1]]] <- paste(kv[-1], collapse = '=') }
req <- function(k) { if (is.null(opt[[k]])) stop('missing --', k); opt[[k]] }
indir <- req('indir'); outdir <- req('outdir'); source(req('rules')); variants <- strsplit(req('variants'), ',', fixed = TRUE)[[1]]; floor_s08 <- as.integer(req('floor'))
if (dir.exists(outdir)) stop('append-only: output dir exists; never overwrite')
dir.create(file.path(outdir, 'tables'), recursive = TRUE)
ages_task <- c(18L, 24L, 30L)

make_design <- function(meta, sex = FALSE) {
    group <- factor(paste(meta$age_weeks, meta$genotype, sep = '_')); x <- model.matrix(~0 + group)
    if (sex && length(unique(meta$sex)) == 2L) x <- cbind(x, sex_M = as.numeric(meta$sex == 'M'))
    if (qr(x)$rank < ncol(x) || nrow(x) <= ncol(x)) stop('unsupported design'); x
}
contrast_vec <- function(x, a) {
    ki <- match(paste0('group', a, '_KI'), colnames(x)); wt <- match(paste0('group', a, '_WT'), colnames(x))
    if (anyNA(c(ki, wt))) return(NULL); v <- numeric(ncol(x)); v[ki] <- 1; v[wt] <- -1; v
}
# Fit one population/age set. Literal TRUE/FALSE only in eBayes (avoids the promise/frame defect found in pipeline_v02).
fit_tables <- function(y, meta, ages, sex = FALSE, trend = FALSE, keep_fit = FALSE) {
    x <- make_design(meta, sex); cms <- lapply(ages, function(a) contrast_vec(x, a)); good <- !vapply(cms, is.null, logical(1))
    if (!any(good)) return(NULL)
    cm <- do.call(cbind, cms[good]); colnames(cm) <- paste0('KI_WT_', ages[good])
    if (trend) { inp <- edgeR::cpm(y, log = TRUE, prior.count = 3); v <- NULL } else { v <- limma::voom(y, x, plot = FALSE); inp <- v }
    fit <- limma::contrasts.fit(limma::lmFit(inp, x), cm)
    fit <- if (trend) limma::eBayes(fit, trend = TRUE, robust = TRUE) else limma::eBayes(fit, robust = TRUE)
    tabs <- lapply(seq_len(ncol(cm)), function(k) limma::topTable(fit, coef = k, number = Inf, confint = 0.95, sort.by = 'none'))
    names(tabs) <- ages[good]
    list(tables = tabs, ages = ages[good], design = x, v = v, fit = fit, cm = cm, df_prior = fit$df.prior, s2_prior = fit$s2.prior)
}

read_counts <- function(file) { r <- read.delim(gzfile(file), check.names = FALSE); ids <- r$stable_id; r$stable_id <- NULL; m <- as.matrix(r); rownames(m) <- ids; m }
meta_all <- read.delim(file.path(indir, 'pb_columns.tsv'), colClasses = c(animal_id = 'character', column = 'character', genotype = 'character', sex = 'character', variant = 'character'))
index_rows <- list(); colrows <- list()
tab_name <- function(v, fl, fit, a) paste0('table__', v, '__floor', fl, '__', fit, '__age', a, '.tsv.gz')
write_tab <- function(tt, v, fl, fit, a, valid) {
    d <- data.frame(stable_id = rownames(tt), logFC = tt$logFC, CI.L = tt$CI.L, CI.R = tt$CI.R, t = tt$t, P.Value = tt$P.Value, AveExpr = tt$AveExpr, stringsAsFactors = FALSE)
    d$BH_within_table <- p.adjust(d$P.Value, 'BH')
    if (!valid) { d$CI.L <- NA_real_; d$CI.R <- NA_real_; d$t <- NA_real_; d$P.Value <- NA_real_; d$BH_within_table <- NA_real_ }   # one KI mouse = estimate only, no p/BH/CI (S08 rule)
    f <- file.path(outdir, 'tables', tab_name(v, fl, fit, a)); stopifnot(!file.exists(f)); con <- gzfile(f, 'wt'); write.table(d, con, sep = '\t', quote = FALSE, row.names = FALSE); close(con)
    basename(f)
}
for (v in variants) {
    mv <- meta_all[meta_all$variant == v, ]; counts <- read_counts(file.path(indir, paste0('pb_', v, '.tsv.gz'))); stopifnot(all(mv$column %in% colnames(counts)))
    for (fl in c(1L, floor_s08)) {
        mm <- mv[mv$n_nuclei_in_column >= fl, , drop = FALSE]; stopifnot(nrow(mm) >= 6L)
        raw <- counts[, mm$column, drop = FALSE]
        y <- prepare_pseudobulk(raw); keep <- expression_filter(y, mm$age_weeks)
        full_sizes <- y$samples$lib.size; nf <- y$samples$norm.factors
        y <- y[keep, , keep.lib.sizes = TRUE]; stopifnot(all(y$samples$lib.size == full_sizes), all(y$samples$norm.factors == nf))
        colrows[[length(colrows) + 1L]] <- data.frame(variant = v, floor = fl, column = mm$column, age_weeks = mm$age_weeks, genotype = mm$genotype, n_nuclei_in_column = mm$n_nuclei_in_column, lib.size = full_sizes, norm.factors = nf)
        mm$age_weeks <- as.integer(mm$age_weeks)
        pooled <- tryCatch(fit_tables(y, mm, ages_task), error = function(e) NULL)
        for (a in ages_task) {
            nw <- sum(mm$age_weeks == a & mm$genotype == 'WT'); nk <- sum(mm$age_weeks == a & mm$genotype == 'KI'); ki_ids <- paste(mm$animal_id[mm$age_weeks == a & mm$genotype == 'KI'], collapse = ',')
            # per-age 3v3 (the S08 decision path for Purkinje): subset with keep.lib.sizes = TRUE
            pa <- NULL; st <- 'NOT_FITTABLE'
            if (nk >= 1L && nw >= 1L) { s <- mm$age_weeks == a; pa <- tryCatch(fit_tables(y[, s, keep.lib.sizes = TRUE], mm[s, , drop = FALSE], a), error = function(e) NULL)
                st <- if (is.null(pa) || is.null(pa$tables[[as.character(a)]])) 'NOT_FITTABLE' else if (min(nw, nk) < 2L) 'ESTIMATE_ONLY' else 'FITTED' }
            tf <- if (st == 'NOT_FITTABLE') NA_character_ else write_tab(pa$tables[[as.character(a)]], v, fl, 'perage', a, st == 'FITTED')
            index_rows[[length(index_rows) + 1L]] <- data.frame(variant = v, floor = fl, fit = 'perage', age = a, status = st, n_WT = nw, n_KI = nk, n_genes = nrow(y), ki_mice = ki_ids, table_file = tf, stringsAsFactors = FALSE)
            stp <- 'NOT_FITTABLE'
            if (!is.null(pooled) && !is.null(pooled$tables[[as.character(a)]])) stp <- if (min(nw, nk) < 2L) 'ESTIMATE_ONLY' else 'FITTED'
            tfp <- if (stp == 'NOT_FITTABLE') NA_character_ else write_tab(pooled$tables[[as.character(a)]], v, fl, 'pooled', a, stp == 'FITTED')
            index_rows[[length(index_rows) + 1L]] <- data.frame(variant = v, floor = fl, fit = 'pooled', age = a, status = stp, n_WT = nw, n_KI = nk, n_genes = nrow(y), ki_mice = ki_ids, table_file = tfp, stringsAsFactors = FALSE)
        }
        cat('fits done:', v, 'floor', fl, '\n'); flush.console()
    }
}
write.table(do.call(rbind, index_rows), file.path(outdir, 'FIT_INDEX.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, colrows), file.path(outdir, 'COLUMNS_TMM.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
writeLines(c(R.version.string, paste('limma', as.character(packageVersion('limma'))), paste('edgeR', as.character(packageVersion('edgeR')))), file.path(outdir, 'VERSIONS.txt'))
writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))
cat('far-set fits complete ->', outdir, '\n')

#%% debug
