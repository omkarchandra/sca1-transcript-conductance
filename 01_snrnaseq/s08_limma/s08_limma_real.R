# S08 real-data limma pipeline (Track A). Ports the frozen edgeR S08 structure (s08_fit_expression.R, s08_robustness.R, s08_boundary_*.R, s08_finalize.R) with limma as the fitter, as ruled in memo v13 (S8-A decision):
#   * primary = voom(frozen effective library sizes) -> lmFit -> contrasts.fit -> eBayes(robust=TRUE) -> topTable(confint=0.95),
#     LITERAL 8-column ~0+genotype_age design (no sex term); per-age KI-WT contrasts by EXACT column-name equality (no regex);
#   * sensitivities: limma-trend (cpm log prior.count=3, eBayes trend+robust); common-sex-term model (where both sexes exist);
#     doublet-excluded refit (frozen doublet_y); 24 leave-one-mouse-out refits; exact within-sex permutation enumeration per age
#     (|moderated t|, descriptive, NO BH); 5W Microglia 1-vs-3 = ESTIMATE_ONLY (no p/BH/CI);
#   * normalisation: the cache's DGEList lib.size x norm.factors ARE the frozen effective sizes; calcNormFactors is NEVER called;
#   * ONE joint BH 0.05 per family (core_primary, whole_transcriptome_secondary) across ALL populations and ages; non-finite p counts
#     as 1 and sets failed_fit_p1; stop if >1% unresolved in any population x age in either family;
#   * order: STAGE 1 (WT-only variance trigger + WT-only MDS, own files + hash manifest) runs BEFORE any KI column enters a fit;
#     STAGE 2 refuses to start unless the manifest exists and matches. Variance trigger is PER POPULATION; fallback (i) per-age
#     fits is computed ONLY for a population that triggers; with --fallback_decides=true (DEFAULT, plain reading of 'engages') the per-age
#     p-values feed the joint BH and the decision/reason of that population's rows, the pooled primary is reported flagged alongside
#     (columns *_pooled / BH_q_pooled); with =false the pooled primary decides and the per-age fit is a side report. An unsupported
#     fallback for a triggered population STOPS the run. This is ruling (e) in the review packet;
#   * camera (limma::camera, inter.gene.cor 0.01) for M1-M10 + Magenta in Purkinje cell (primary) and Bergmann glia (secondary),
#     per age contrast, BH across modules within each contrast; <5 tested genes = NOT_TESTABLE; absent genes listed;
#   * output = one row per registry row (27,415 core_primary + 531,565 secondary), S11 vocabulary: stable_id, population, age
#     (integer; age_weeks kept = schema delta), decision, reason, failed_fit_p1.
# Inputs are ONLY paths + sha256 in --pins (JSON); any mismatch stops. --mode=real needs --approval_file naming this script's sha256.
# Usage: Rscript s08_limma_real_v01_A.R --out=DIR --pins=pins.json --mode=synthetic|real [--approval_file=F] [--n_perm_variance=10000]
#        [--seed=20261006] [--fallback_decides=false] [--stage=all|gates] [--force_trigger_synthetic=POP  (synthetic mode only)]
suppressPackageStartupMessages({library(limma); library(edgeR); library(jsonlite)})
`%||%` <- function(a, b) if (is.null(a)) b else a
args <- commandArgs(trailingOnly = TRUE); opt <- list()
for (a in args) { kv <- strsplit(sub('^--', '', a), '=', fixed = TRUE)[[1]]; opt[[kv[1]]] <- paste(kv[-1], collapse = '=') }
req <- function(k) { if (is.null(opt[[k]])) stop('missing --', k); opt[[k]] }
script_path <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1])
sha <- function(p) strsplit(system2('sha256sum', shQuote(p), stdout = TRUE), ' ')[[1]][1]
out <- req('out'); mode <- req('mode'); stopifnot(mode %in% c('synthetic', 'real'))
seed <- as.integer(opt[['seed']] %||% 20261006L); B_var <- as.integer(opt[['n_perm_variance']] %||% 10000L)
fallback_decides <- !identical(tolower(opt[['fallback_decides']] %||% 'true'), 'false'); stage <- opt[['stage']] %||% 'all'
s11_family <- opt[['s11_family']] %||% 'core_primary'; t_start <- proc.time()[['elapsed']]; timings <- list()
if (dir.exists(out)) stop('append-only: output dir exists; never overwrite')
pins <- read_json(req('pins'))
pin <- function(name) { p <- pins[[name]]; if (is.null(p) || !file.exists(p$path)) stop('pin missing: ', name); if (!identical(sha(p$path), p$sha256)) stop('sha256 mismatch for pin: ', name); p$path }
paths <- sapply(c('cache', 'registry', 'lock', 'prefit', 'ambient', 'contrasts', 'dictionary', 'module_ids', 'memo', 'decision_record'), pin)
approval <- NULL
if (mode == 'real') {
    ap <- req('approval_file'); stopifnot(file.exists(ap)); approval <- read_json(ap)
    stopifnot(identical(approval$cleared_script_sha256, sha(script_path)), identical(approval$memo_sha256, pins$memo$sha256), identical(approval$decision_record_sha256, pins$decision_record$sha256))
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)
lock <- read_json(paths[['lock']], simplifyVector = TRUE)
stopifnot(!lock$biological_fits_run %||% FALSE || mode == 'synthetic')
registry <- read.delim(gzfile(paths[['registry']]), stringsAsFactors = FALSE)
cache <- readRDS(paths[['cache']])
ambient <- read.delim(paths[['ambient']], stringsAsFactors = FALSE)
core_meta <- read.delim(paths[['dictionary']], stringsAsFactors = FALSE)
module_ids <- read_json(paths[['module_ids']], simplifyVector = TRUE)$modules
stopifnot(all(c('stable_id', 'population', 'age_weeks', 'family') %in% names(registry)), !anyDuplicated(registry[, c('stable_id', 'population', 'age_weeks', 'family')]))
if (mode == 'real') stopifnot(sum(registry$family == 'core_primary') == 27415L, sum(registry$family != 'core_primary') == 531565L)
pops <- names(cache); reg_pops <- sort(unique(registry$population)); stopifnot(all(reg_pops %in% pops))

# ------------------------------------------------------------------ STAGE 1: WT-only gates (no KI column is used here)
stage1_dir <- file.path(out, 'stage1_wt_only'); dir.create(stage1_dir)
gate_rows <- list(); mds_rows <- list()
variance_gate <- function(pop, object) {
    y <- object$primary_y; meta <- object$primary_meta
    wt <- which(meta$genotype == 'WT'); stopifnot(length(wt) >= 6L)           # KI columns are never subset below this line
    yw <- y[, wt, keep.lib.sizes = TRUE]; age <- meta$age_weeks[wt]
    lc <- edgeR::cpm(yw, log = TRUE, prior.count = 3)                          # normalised (frozen effective) library sizes
    ages <- sort(unique(age)); keep_age <- ages[table(factor(age, levels = ages)) >= 2L]; stopifnot(length(keep_age) >= 2L)
    sel <- age %in% keep_age; lc <- lc[, sel, drop = FALSE]; age <- age[sel]
    stat_T <- function(lab, genes = NULL) {
        M <- model.matrix(~0 + factor(lab)); na <- colSums(M); G <- nrow(lc)
        S1 <- lc %*% M; S2 <- (lc^2) %*% M
        v <- (S2 - S1^2 / matrix(na, G, ncol(M), byrow = TRUE)) / matrix(na - 1, G, ncol(M), byrow = TRUE)
        lv <- log(pmax(v, 1e-10)); A <- ncol(lv); sv <- (rowSums(lv^2) - rowSums(lv)^2 / A) / (A - 1)
        list(T = median(if (is.null(genes)) sv else sv[genes]), sv = sv, v = v)
    }
    obs <- stat_T(age); ok <- rowSums(obs$v <= 1e-10) == 0; Tobs <- median(obs$sv[ok])
    set.seed(seed)
    Tnull <- vapply(seq_len(B_var), function(b) stat_T(sample(age), ok)$T, numeric(1))
    q <- quantile(Tnull, c(.025, .5, .95, .975), names = FALSE)
    list(row = data.frame(population = pop, n_WT = length(wt), ages_used = paste(keep_age, collapse = ','), n_genes_used = sum(ok), n_genes_excluded_zero_variance = sum(!ok),
        T_observed = Tobs, null_q025 = q[1], null_q50 = q[2], null_q95 = q[3], null_q975 = q[4], p_perm = (1 + sum(Tnull >= Tobs)) / (B_var + 1L),
        B = B_var, seed = seed, triggered = Tobs > q[3], stringsAsFactors = FALSE),
        mds = { m <- limma::plotMDS(lc, plot = FALSE); data.frame(population = pop, animal_id = meta$animal_id[wt][sel], age_weeks = age, sex = meta$sex[wt][sel], dim1 = m$x, dim2 = m$y) })
}
for (pop in reg_pops) { g <- variance_gate(pop, cache[[pop]]); gate_rows[[pop]] <- g$row; mds_rows[[pop]] <- g$mds }
gates <- do.call(rbind, gate_rows); mds <- do.call(rbind, mds_rows)
gates$trigger_forced_synthetic <- FALSE
if (mode == 'synthetic' && !is.null(opt[['force_trigger_synthetic']])) { f <- gates$population == opt[['force_trigger_synthetic']]; gates$triggered[f] <- TRUE; gates$trigger_forced_synthetic[f] <- TRUE }   # synthetic-only switch to exercise fallback (i); ignored in real mode
write.table(gates, file.path(stage1_dir, 'WT_ONLY_VARIANCE_TRIGGER.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(mds, file.path(stage1_dir, 'WT_ONLY_MDS.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
manifest <- list(stage = 'wt_only_gates_complete', utc = format(Sys.time(), tz = 'UTC', usetz = TRUE), script_sha256 = sha(script_path),
    files = list(WT_ONLY_VARIANCE_TRIGGER.tsv = sha(file.path(stage1_dir, 'WT_ONLY_VARIANCE_TRIGGER.tsv')), WT_ONLY_MDS.tsv = sha(file.path(stage1_dir, 'WT_ONLY_MDS.tsv'))),
    ki_columns_used = FALSE, triggered_populations = gates$population[gates$triggered])
write_json(manifest, file.path(stage1_dir, 'STAGE1_MANIFEST.json'), auto_unbox = TRUE, pretty = TRUE)
timings$stage1_seconds <- proc.time()[['elapsed']] - t_start
cat('stage 1 complete; triggered populations:', if (any(gates$triggered)) paste(gates$population[gates$triggered], collapse = ', ') else 'none', '\n')
if (stage == 'gates') quit(status = 0)

# ------------------------------------------------------------------ STAGE 2: refuses without a matching stage-1 manifest
mf <- read_json(file.path(stage1_dir, 'STAGE1_MANIFEST.json'))
stopifnot(identical(mf$stage, 'wt_only_gates_complete'), identical(unlist(mf$files), c(WT_ONLY_VARIANCE_TRIGGER.tsv = sha(file.path(stage1_dir, 'WT_ONLY_VARIANCE_TRIGGER.tsv')), WT_ONLY_MDS.tsv = sha(file.path(stage1_dir, 'WT_ONLY_MDS.tsv')))))
triggered <- gates$population[gates$triggered]

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
boundary_type <- function(wt, ki) ifelse(wt == 0 & ki == 0, 'not_identified_both_genotypes_zero', ifelse(wt == 0, 'positive_boundary', ifelse(ki == 0, 'negative_boundary', 'finite')))
fam_levels <- unique(registry$family)

res <- registry; res$age <- as.integer(res$age_weeks)
for (cname in c('log2FC', 'CI_lower', 'CI_upper', 't', 'P_value', 'AveExpr', 'WT_raw_sum', 'KI_raw_sum', 'trend_log2FC', 'trend_P', 'sex_log2FC', 'sex_P', 'doublet_log2FC', 'doublet_P', 'fallback_perage_log2FC', 'fallback_perage_P', 'fallback_perage_CI_lower', 'fallback_perage_CI_upper')) res[[cname]] <- NA_real_
res$n_WT <- 0L; res$n_KI <- 0L; res$failure <- ''; res$variance_trigger_flag <- res$population %in% triggered
estimate_only <- list(); diagnostics <- list(); fits_primary <- list(); camera_rows <- list(); fit_errors <- list()
try_fit <- function(label, expr) tryCatch(expr, error = function(e) { fit_errors[[label]] <<- conditionMessage(e); NULL })
for (pop in reg_pops) {
    obj <- cache[[pop]]; y <- obj$primary_y; mm <- obj$primary_meta; core <- obj$core_idx
    planned <- which(res$population == pop); ages <- sort(unique(res$age_weeks[planned]))
    prim <- try_fit(paste(pop, 'primary'), fit_tables(y, mm, ages))
    trd <- try_fit(paste(pop, 'trend'), fit_tables(y, mm, ages, trend = TRUE))
    sx <- try_fit(paste(pop, 'sex_term'), fit_tables(y, mm, ages, sex = TRUE))
    dbl <- try_fit(paste(pop, 'doublet'), fit_tables(obj$doublet_y, obj$doublet_meta, ages))
    fb <- if (pop %in% triggered) lapply(ages, function(a) try_fit(paste(pop, 'fallback_age', a), { s <- mm$age_weeks == a; fit_tables(y[, s, keep.lib.sizes = TRUE], mm[s, , drop = FALSE], a) })) else NULL
    if (!is.null(fb) && (any(vapply(fb, is.null, logical(1))) || any(vapply(seq_along(ages), function(i) is.null(fb[[i]]$tables[[as.character(ages[i])]]), logical(1))))) stop('fallback (i) per-age fits cannot be supported for triggered population ', pop, ': STOP per ruling (report, do not switch methods)')
    fits_primary[[pop]] <- list(prim = prim)
    for (a in ages) {
        ii <- planned[res$age_weeks[planned] == a]; ak <- as.character(a)
        res$n_WT[ii] <- sum(mm$age_weeks == a & mm$genotype == 'WT'); res$n_KI[ii] <- sum(mm$age_weeks == a & mm$genotype == 'KI')
        if (is.null(prim) || is.null(prim$tables[[ak]])) { res$failure[ii] <- 'population_fit_failed_or_contrast_missing'; next }
        tt <- prim$tables[[ak]]; ix <- match(res$stable_id[ii], rownames(tt)); stopifnot(!anyNA(ix))
        res$log2FC[ii] <- tt$logFC[ix]; res$CI_lower[ii] <- tt$CI.L[ix]; res$CI_upper[ii] <- tt$CI.R[ix]; res$t[ii] <- tt$t[ix]
        res$P_value[ii] <- tt$P.Value[ix]; res$AveExpr[ii] <- tt$AveExpr[ix]
        wt_sum <- rowSums(y$counts[, mm$age_weeks == a & mm$genotype == 'WT', drop = FALSE]); ki_sum <- rowSums(y$counts[, mm$age_weeks == a & mm$genotype == 'KI', drop = FALSE])
        res$WT_raw_sum[ii] <- wt_sum[ix]; res$KI_raw_sum[ii] <- ki_sum[ix]
        for (item in list(list('trend', trd, 'trend_log2FC', 'trend_P'), list('sex', sx, 'sex_log2FC', 'sex_P'), list('doublet', dbl, 'doublet_log2FC', 'doublet_P'))) {
            f <- item[[2]]; if (is.null(f) || is.null(f$tables[[ak]])) next
            if (item[[1]] == 'doublet' && min(sum(obj$doublet_meta$age_weeks == a & obj$doublet_meta$genotype == 'WT'), sum(obj$doublet_meta$age_weeks == a & obj$doublet_meta$genotype == 'KI')) < 2L) next
            t2 <- f$tables[[ak]]; j <- match(res$stable_id[ii], rownames(t2)); ok <- !is.na(j)
            res[[item[[3]]]][ii[ok]] <- t2$logFC[j[ok]]; res[[item[[4]]]][ii[ok]] <- t2$P.Value[j[ok]]
        }
        if (!is.null(fb)) { f <- fb[[match(a, ages)]]; if (!is.null(f) && !is.null(f$tables[[ak]])) { t3 <- f$tables[[ak]]; j <- match(res$stable_id[ii], rownames(t3)); ok <- !is.na(j)
            res$fallback_perage_log2FC[ii[ok]] <- t3$logFC[j[ok]]; res$fallback_perage_P[ii[ok]] <- t3$P.Value[j[ok]]; res$fallback_perage_CI_lower[ii[ok]] <- t3$CI.L[j[ok]]; res$fallback_perage_CI_upper[ii[ok]] <- t3$CI.R[j[ok]] } }
    }
    if (!is.null(prim)) diagnostics[[pop]] <- list(n_mice = ncol(y), n_genes = nrow(y), design_rank = qr(prim$design)$rank, residual_df = nrow(mm) - ncol(prim$design),
        df_prior_median = median(prim$df_prior), df_prior_min = min(prim$df_prior), s2_prior_median = median(prim$s2_prior), variance_trigger = pop %in% triggered)
    if (!is.null(prim)) {                                                   # an age present in the design but not planned, with 1 mouse in a genotype = ESTIMATE_ONLY (no p/BH/CI)
        for (a in setdiff(sort(unique(mm$age_weeks)), ages)) {
            nw <- sum(mm$age_weeks == a & mm$genotype == 'WT'); nk <- sum(mm$age_weeks == a & mm$genotype == 'KI'); cv <- contrast_vec(prim$design, a)
            if (min(nw, nk) >= 1L && !is.null(cv)) { ff <- limma::lmFit(prim$v, prim$design); eff <- drop(ff$coefficients %*% cv)
                estimate_only[[paste(pop, a)]] <- data.frame(stable_id = rownames(y)[core], population = pop, age = as.integer(a), age_weeks = a, log2FC = eff[core], n_WT = nw, n_KI = nk,
                    status = 'ESTIMATE_ONLY', reason = 'one_vs_n_estimate_only_no_p_BH_CI', stringsAsFactors = FALSE) }
        }
    }
    cat('stage 2 fits complete:', pop, '\n'); flush.console()
}
# ---- decision-driving fields: pooled primary, or the per-age fallback for triggered populations (switch)
res$decision_source <- 'pooled_primary'; res$P_decision <- res$P_value; res$log2FC_decision <- res$log2FC; res$CI_lower_decision <- res$CI_lower; res$CI_upper_decision <- res$CI_upper
if (fallback_decides && length(triggered)) { f <- res$population %in% triggered; res$decision_source[f] <- 'fallback_per_age'; res$P_decision[f] <- res$fallback_perage_P[f]
    res$log2FC_decision[f] <- res$fallback_perage_log2FC[f]; res$CI_lower_decision[f] <- res$fallback_perage_CI_lower[f]; res$CI_upper_decision[f] <- res$fallback_perage_CI_upper[f] }
# ---- BH per family over ALL populations and ages; non-finite p = 1 and failed_fit_p1
res$failed_fit_p1 <- !is.finite(res$P_decision); res$BH_q <- NA_real_; res$BH_q_pooled <- NA_real_
for (fam in fam_levels) { ii <- which(res$family == fam); res$BH_q[ii] <- p.adjust(ifelse(is.finite(res$P_decision[ii]), res$P_decision[ii], 1), 'BH')
    res$BH_q_pooled[ii] <- p.adjust(ifelse(is.finite(res$P_value[ii]), res$P_value[ii], 1), 'BH') }
for (pre in c('trend', 'sex', 'doublet')) {
    p <- res[[paste0(pre, '_P')]]; q <- rep(NA_real_, nrow(res))
    for (fam in fam_levels) { ii <- which(res$family == fam); q[ii] <- p.adjust(ifelse(is.finite(p[ii]), p[ii], 1), 'BH') }
    res[[paste0(pre, '_BH_q')]] <- q; l <- res[[paste0(pre, '_log2FC')]]
    res[[paste0(pre, '_opposite_sign')]] <- is.finite(l) & is.finite(res$log2FC) & l * res$log2FC < 0
    res[[paste0(pre, '_BH_call_change')]] <- (q < .05) != (res$BH_q < .05)
    res[[paste0(pre, '_eligibility_change')]] <- is.finite(p) != is.finite(res$P_value)
}
res$effect_status <- boundary_type(res$WT_raw_sum, res$KI_raw_sum)
akey <- paste(ambient$stable_id, ambient$target_population, sep = '::'); res$potential_ambient_suspect <- paste(res$stable_id, res$population, sep = '::') %in% akey
res$ambient_flagged_hit <- res$potential_ambient_suspect & res$BH_q < .05
res$core_canonical_symbols <- core_meta$core_symbols[match(res$stable_id, core_meta$stable_id)]
res$decision <- ifelse(res$failed_fit_p1, 'UNRESOLVED_OR_WITHHELD', ifelse(res$BH_q < .05, 'HIT', 'TESTED_NOT_A_HIT'))
res$reason <- paste0(ifelse(res$decision_source == 'fallback_per_age', 'fallback_per_age:', ''), ifelse(res$failed_fit_p1, 'non_finite_fit_counted_as_p1', ifelse(res$BH_q < .05, 'family_BH_lt_0.05', 'valid_test_BH_ge_0.05')))
# ---- stop rule per family x population x age
unres <- aggregate(failed_fit_p1 ~ family + population + age, data = res, FUN = mean); unres$stop <- unres$failed_fit_p1 > 0.01
write.table(unres, file.path(out, 'stop_rule_unresolved_fractions.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
if (any(unres$stop)) { write_json(list(stop = TRUE, cells = unres[unres$stop, ]), file.path(out, 'STOP_RULE.json'), auto_unbox = TRUE, pretty = TRUE)
    write.table(res, gzfile(file.path(out, 'expression_results_PARTIAL_STOPPED.tsv.gz')), sep = '\t', quote = FALSE, row.names = FALSE); stop('stop rule: >1% unresolved in some family x population x age; see STOP_RULE.json') }

timings$stage2_seconds <- proc.time()[['elapsed']] - t_start
# ---- STAGE 3: sensitivities on the core family (LOMO refits; exact within-sex permutation enumeration)
core_reg <- res[res$family == 'core_primary', ]
cartesian_assignments <- function(m, a) {
    rows <- which(m$age_weeks == a); strata <- split(rows, m$sex[rows])
    choices <- lapply(strata, function(z) { k <- sum(m$genotype[z] == 'KI'); if (k == 0) return(list(integer(0))); if (k == length(z)) return(list(z)); combn(z, k, simplify = FALSE) })
    grid <- expand.grid(lapply(choices, seq_along))
    lapply(seq_len(nrow(grid)), function(i) { g <- m$genotype; g[rows] <- 'WT'; for (j in seq_along(choices)) g[choices[[j]][[grid[i, j]]]] <- 'KI'; g })
}
perm_list <- list(); lomo_list <- list(); alloc_list <- list()
for (pop in reg_pops) {
    obj <- cache[[pop]]; y <- obj$primary_y; m <- obj$primary_meta; core <- obj$core_idx
    ages <- sort(unique(core_reg$age_weeks[core_reg$population == pop])); prim <- fits_primary[[pop]]$prim
    for (a in ages) {
        ak <- as.character(a); obs_t <- prim$tables[[ak]]$t[core]
        asg <- cartesian_assignments(m, a); numer <- numeric(length(core)); failed <- 0L
        for (g in asg) {
            mm <- m; mm$genotype <- g
            tv <- tryCatch({ f <- fit_tables(y, mm, a); f$tables[[ak]]$t[core] }, error = function(e) rep(NA_real_, length(core)))
            if (any(!is.finite(tv))) failed <- failed + 1L
            numer <- numer + ifelse(is.finite(tv), abs(tv) >= abs(obs_t) - 1e-10 * pmax(1, abs(obs_t)), 0)
        }
        Bn <- length(asg)
        perm_list[[length(perm_list) + 1L]] <- data.frame(stable_id = rownames(y)[core], population = pop, age_weeks = a, numerator = numer, B = Bn,
            p_exact = if (failed == 0L && Bn > 1L) numer / Bn else NA_real_, failed_assignments = failed, generic_bound = 1 / Bn,
            status = if (Bn > 1L) 'descriptive_secondary_robustness_no_BH' else 'unavailable_no_informative_allocation', stringsAsFactors = FALSE)
        alloc_list[[length(alloc_list) + 1L]] <- data.frame(population = pop, age_weeks = a, B = Bn, failed_assignments = failed)
    }
    for (d in seq_len(nrow(m))) {
        mm <- m[-d, , drop = FALSE]; yy <- y[, -d, keep.lib.sizes = TRUE]
        f <- tryCatch(fit_tables(yy, mm, ages), error = function(e) NULL)
        for (a in ages) {
            nw <- sum(mm$age_weeks == a & mm$genotype == 'WT'); nk <- sum(mm$age_weeks == a & mm$genotype == 'KI'); ak <- as.character(a)
            val <- if (is.null(f) || is.null(f$tables[[ak]]) || min(nw, nk) == 0L) rep(NA_real_, length(core)) else f$tables[[ak]]$logFC[match(rownames(y)[core], rownames(f$tables[[ak]]))]
            lomo_list[[length(lomo_list) + 1L]] <- data.frame(stable_id = rownames(y)[core], population = pop, age_weeks = a, dropped_mouse = as.character(m$animal_id[d]), log2FC = val,
                n_WT = nw, n_KI = nk, status = if (all(is.finite(val))) { if (min(nw, nk) == 1L) 'singleton_sensitivity_only' else 'supported_sensitivity' } else 'unsupported', stringsAsFactors = FALSE)
        }
    }
    cat('stage 3 sensitivities complete:', pop, '\n'); flush.console()
}
perm <- do.call(rbind, perm_list); lomo <- do.call(rbind, lomo_list); alloc <- do.call(rbind, alloc_list)
key <- function(d) paste(d$stable_id, d$population, d$age_weeks, sep = '::')
fullp <- res[res$family == 'core_primary', ]; full <- fullp$log2FC[match(key(lomo), key(fullp))]
lomo$opposite_sign <- is.finite(full) & is.finite(lomo$log2FC) & sign(full) * sign(lomo$log2FC) < 0
lomo$difference_from_full <- ifelse(is.finite(full) & is.finite(lomo$log2FC), lomo$log2FC - full, NA_real_)
grp <- split(seq_len(nrow(lomo)), key(lomo))
ls <- do.call(rbind, lapply(names(grp), function(k) { d <- lomo[grp[[k]], ]; fin <- d$difference_from_full[is.finite(d$difference_from_full)]
    data.frame(key = k, LOMO_opposite_sign = any(d$opposite_sign), LOMO_unsupported = sum(d$status == 'unsupported'), LOMO_min_difference = if (length(fin)) min(fin) else NA_real_, LOMO_max_difference = if (length(fin)) max(fin) else NA_real_) }))
li <- match(key(res), ls$key); for (cn in setdiff(names(ls), 'key')) res[[cn]] <- ls[[cn]][li]
pidx <- match(key(res), key(perm)); for (cn in c('numerator', 'B', 'p_exact', 'failed_assignments', 'generic_bound')) res[[paste0('permutation_', cn)]] <- perm[[cn]][pidx]

timings$stage3_seconds <- proc.time()[['elapsed']] - t_start
# ---- STAGE 4: camera for modules (Purkinje cell primary, Bergmann glia secondary)
for (pop in intersect(c('Purkinje cell', 'Bergmann glia'), reg_pops)) {
    obj <- cache[[pop]]; y <- obj$primary_y; prim <- fits_primary[[pop]]$prim; if (is.null(prim)) next
    for (mod in names(module_ids)) {
        spec <- module_ids[[mod]]; if (!(pop %in% spec$populations)) next
        idx <- which(rownames(y) %in% spec$stable_ids); absent <- setdiff(spec$stable_ids, rownames(y))
        for (a in prim$ages) {
            k <- match(a, prim$ages)
            if (length(idx) < 5L) { camera_rows[[length(camera_rows) + 1L]] <- data.frame(population = pop, module = mod, age = as.integer(a), n_tested = length(idx), n_absent = length(absent), Direction = NA, PValue = NA_real_, status = 'NOT_TESTABLE_lt5', stringsAsFactors = FALSE); next }
            cr <- tryCatch(limma::camera(prim$v, index = idx, design = prim$design, contrast = prim$cm[, k], inter.gene.cor = 0.01), error = function(e) NULL)
            camera_rows[[length(camera_rows) + 1L]] <- data.frame(population = pop, module = mod, age = as.integer(a), n_tested = length(idx), n_absent = length(absent),
                Direction = if (is.null(cr)) NA else as.character(cr$Direction), PValue = if (is.null(cr)) NA_real_ else cr$PValue, status = if (is.null(cr)) 'camera_failed' else 'TESTED', stringsAsFactors = FALSE)
        }
    }
}
cam <- do.call(rbind, camera_rows); cam$BH_q_within_contrast <- NA_real_
for (g in split(seq_len(nrow(cam)), paste(cam$population, cam$age))) { p <- cam$PValue[g]; ok <- is.finite(p); cam$BH_q_within_contrast[g[ok]] <- p.adjust(p[ok], 'BH') }
cam$expected_direction <- vapply(cam$module, function(m) module_ids[[m]]$expected_direction, character(1))
cam$expected_direction_label <- ifelse(cam$status != 'TESTED', NA, ifelse(cam$expected_direction == 'none', 'no_expectation', ifelse(tolower(cam$Direction) == cam$expected_direction, 'matches', 'contradicts')))
cam$role <- ifelse(cam$population == 'Purkinje cell', 'primary', 'secondary')

timings$stage4_seconds <- proc.time()[['elapsed']] - t_start
# ---- outputs
write.table(res, gzfile(file.path(out, 'expression_results_complete.tsv.gz')), sep = '\t', quote = FALSE, row.names = FALSE)
# S11 release table: S11 (s11_fragments_v01.gene_node) accepts only populations 'purkinje' / 'bergmann' and ages 5/18/24/30 -> release subset with those codes
s11 <- res[res$population %in% c('Purkinje cell', 'Bergmann glia') & res$family == s11_family, ]
stopifnot(!anyDuplicated(s11[, c('stable_id', 'population', 'age')])); s11$population_registry_name <- s11$population; s11$population <- ifelse(s11$population == 'Purkinje cell', 'purkinje', 'bergmann')
write.table(s11, gzfile(file.path(out, 'S08_S11_release_purkinje_bergmann.tsv.gz')), sep = '\t', quote = FALSE, row.names = FALSE)
write_json(list(schema_delta = list(added_columns = c('age (integer = age_weeks)', 'decision', 'reason', 'failed_fit_p1', 'effect_status', 'population_registry_name (release table only)'),
    release_table = paste0('S08_S11_release_purkinje_bergmann.tsv.gz holds ONLY family=', s11_family, ' rows (core genes also appear in the secondary family; one row per stable_id x population x age) for Purkinje cell and Bergmann glia with S11 population codes purkinje/bergmann; the complete table keeps registry population names for all 13 populations and both families'),
    decision_vocabulary = c('HIT', 'TESTED_NOT_A_HIT', 'UNRESOLVED_OR_WITHHELD'), estimate_only_file = 'estimate_only.tsv (status ESTIMATE_ONLY)')), file.path(out, 'SCHEMA_DELTA.json'), auto_unbox = TRUE, pretty = TRUE)
write.table(perm, gzfile(file.path(out, 'permutation.tsv.gz')), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(lomo, gzfile(file.path(out, 'LOMO.tsv.gz')), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(alloc, file.path(out, 'allocation_spaces.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(cam, file.path(out, 'module_camera_results.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
if (length(estimate_only)) write.table(do.call(rbind, estimate_only), file.path(out, 'estimate_only.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
diagnostics$fit_errors <- fit_errors
write_json(diagnostics, file.path(out, 'model_diagnostics.json'), auto_unbox = TRUE, pretty = TRUE)
fam_sum <- lapply(fam_levels, function(f) { r <- res[res$family == f, ]; list(family = f, planned = nrow(r), successful = sum(is.finite(r$P_value)), HIT = sum(r$decision == 'HIT'), unresolved = sum(r$failed_fit_p1)) })
write_json(list(stage = 'S08_limma_real_complete', mode = mode, families = fam_sum, variance_triggered_populations = triggered, design = 'literal 8-column ~0+genotype_age per population (no sex term); sex-term model is a sensitivity',
    cross_age_tests = 'not run, descriptive only', one_vs_three = 'estimate only; no p/BH/CI', interpretation = 'none in this file'), file.path(out, 'RESULT.json'), auto_unbox = TRUE, pretty = TRUE)
write_json(list(utc = format(Sys.time(), tz = 'UTC', usetz = TRUE), script_sha256 = sha(script_path), mode = mode, seed = seed, n_perm_variance = B_var,
    pins = lapply(pins, function(p) list(path = p$path, sha256 = p$sha256)), approval_file_sha256 = if (is.null(approval)) 'n/a' else sha(opt[['approval_file']]),
    approval = approval, timings_seconds = c(timings, total = proc.time()[['elapsed']] - t_start), R = R.version.string, limma = as.character(packageVersion('limma')), edgeR = as.character(packageVersion('edgeR')), fallback_decides = fallback_decides),
    file.path(out, 'RUN_RECORD.json'), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
cat('S08 limma real pipeline complete ->', out, '\n')

#%% debug
