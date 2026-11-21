# v03: --root DIR (the A1 package folder) as the first two arguments; the S08 rules are read from DIR/s08_expression_rules.R
# (a copy of ${REPO_ROOT}/01_snrnaseq/s08_limma/s08_expression_rules.R, same hash check); no other change.
# A1 step (c): limma-voom on the Purkinje pseudobulk of each dataset,
# mutant vs control, as S08: TMM on all genes, the S08 expression_filter (sourced from ${REPO_ROOT}/01_snrnaseq/s08_limma/s08_expression_rules.R), voom,
# lmFit on ~0 + genotype (and age where there are several), eBayes(robust = TRUE), topTable with a 95 % confint. The nine genes
# and P2's rule: >= 4 of 5 channel genes with log2FC < 0 and the one-sided sign test over the nine with p < 0.05.
# Reads the NEW real_v01/pseudobulk/ (a1_pseudobulk.py); writes a NEW real_v01/tests/.
#   Rscript ${REPO_ROOT}/02_replication/a1_purkinje/a1_limma.R.
suppressPackageStartupMessages({library(edgeR); library(limma)})
args <- commandArgs(trailingOnly = TRUE)
if (args[1] != '--root') stop('usage: a1_limma.R --root DIR PSEUDOBULK_SUMS [SCA7_SUMS]')
here <- normalizePath(args[2]); args <- args[-(1:2)]
pbdir <- file.path(here, 'real_v01', 'pseudobulk'); out <- file.path(here, 'real_v01', 'tests')
sha <- function(p) strsplit(system2('sha256sum', shQuote(p), stdout = TRUE), ' ')[[1]][1]
rules <- 's08_expression_rules.R'
if (!identical(sha(file.path(here, rules)), 'e1b64c9f6259812cc9d758cea6867f3ef8f4150e1bc2553a46e33d5660449981')) stop('S08 rules hash differs')
source(file.path(here, rules))
pb_sums <- args[1]
if (!identical(sha(file.path(pbdir, 'SHA256SUMS')), pb_sums)) stop('pseudobulk SHA256SUMS differs from the given hash')
if (dir.exists(out)) stop('tests/ already exists')
FIVE <- c('Kcnma1', 'Kcnc3', 'Kcnd3', 'Cacna1a', 'Cacna1g'); NINE <- c(FIVE, 'Itpr1', 'Atp2a3', 'Kcnmb4', 'Cacna2d2')

run <- function(ds, upper) {
  d0 <- if (ds == 'GSE269430_SCA7') sca7 else file.path(pbdir, ds)
  cnt <- as.matrix(read.delim(file.path(d0, 'COUNTS.tsv'), row.names = 1, check.names = FALSE))
  meta <- read.delim(file.path(d0, 'SAMPLES.tsv'), row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if ('kept' %in% names(meta)) meta <- meta[meta$kept %in% c(TRUE, 'True'), ]
  meta <- meta[colnames(cnt), ]
  stopifnot(all(cnt == round(cnt)), nrow(meta) == ncol(cnt))
  y <- calcNormFactors(DGEList(counts = cnt), method = 'TMM')
  keep <- expression_filter(y, meta$age_weeks)
  y <- y[keep, , keep.lib.sizes = TRUE]
  ages <- sort(unique(meta$age_weeks)); res <- list(); reads <- list()
  grp <- factor(if (length(ages) > 1) paste(meta$genotype, meta$age_weeks, sep = '_') else meta$genotype)
  design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
  v <- voom(y, design, plot = FALSE); fit <- lmFit(v, design)
  for (a in ages) {
    lab <- if (length(ages) > 1) c(paste0('KI_', a), paste0('WT_', a)) else c('KI', 'WT')
    cm <- makeContrasts(contrasts = paste(lab[1], '-', lab[2]), levels = design)
    f2 <- eBayes(contrasts.fit(fit, cm), robust = TRUE)
    tt <- topTable(f2, coef = 1, number = Inf, confint = 0.95, sort.by = 'none')
    genes <- if (upper) toupper(NINE) else NINE
    r <- tt[match(genes, rownames(tt)), c('logFC', 'CI.L', 'CI.R', 'P.Value', 'adj.P.Val')]
    r <- data.frame(dataset = ds, age = a, gene = NINE, row_symbol = genes, tested = !is.na(r$logFC), r, row.names = NULL)
    lf <- r$logFC; five <- lf[r$gene %in% FIVE & !is.na(lf)]; nine <- lf[!is.na(lf)]
    sp <- if (length(nine)) binom.test(sum(nine < 0), length(nine), 0.5, alternative = 'greater')$p.value else NA
    reads[[length(reads) + 1]] <- data.frame(dataset = ds, age = a, n_samples_KI = sum(meta$genotype == 'KI' & meta$age_weeks == a),
      n_samples_WT = sum(meta$genotype == 'WT' & meta$age_weeks == a), n_genes_tested = nrow(y), n_channel_found = length(five),
      n_channel_negative = sum(five < 0), n_nine_found = length(nine), n_nine_negative = sum(nine < 0), sign_test_p = sp,
      replicates = sum(five < 0) >= 4 && !is.na(sp) && sp < 0.05, median_channel_log2fc = if (length(five)) median(five) else NA,
      not_tested = paste(r$gene[!r$tested], collapse = ','))
    res[[length(res) + 1]] <- r
  }
  list(nine = do.call(rbind, res), reading = do.call(rbind, reads))
}

datasets <- list(GSE215336 = FALSE, GSE246183_human = TRUE, GSE246183_human_all20 = TRUE)
sca7 <- file.path(here, 'real_v01', 'pseudobulk_sca7')
if (length(args) > 1) {
  if (!identical(sha(file.path(sca7, 'SHA256SUMS')), args[2])) stop('SCA7 pseudobulk SHA256SUMS differs')
  datasets$GSE269430_SCA7 <- FALSE
}
all <- lapply(names(datasets), function(d) run(d, datasets[[d]]))
nine <- do.call(rbind, lapply(all, `[[`, 'nine')); reading <- do.call(rbind, lapply(all, `[[`, 'reading'))
dir.create(out)
for (nm in c('NINE_GENES', 'READING')) {
  con <- file(file.path(out, paste0(nm, '.tsv')), open = 'wx')
  writeLines('# A1: limma-voom on the Purkinje pseudobulk, mutant vs control, S08 settings; P2 rule (>= 4 of 5 negative; one-sided sign test over nine)', con)
  write.table(if (nm == 'NINE_GENES') nine else reading, con, sep = '\t', quote = FALSE, row.names = FALSE); close(con)
}
rec <- list(script_sha256 = sha(file.path(here, 'a1_limma.R')), pseudobulk_sums = pb_sums, s08_rules = rules,
            edgeR = as.character(packageVersion('edgeR')), limma = as.character(packageVersion('limma')), R = R.version.string,
            datasets = names(datasets))
con <- file(file.path(out, 'RECORD.json'), open = 'wx'); writeLines(jsonlite::toJSON(rec, auto_unbox = TRUE, pretty = TRUE, digits = NA), con); close(con)
con <- file(file.path(out, 'SHA256SUMS'), open = 'wx')
for (f in sort(list.files(out))) if (f != 'SHA256SUMS') writeLines(paste0(sha(file.path(out, f)), '  ', f), con)
close(con)
print(reading); print(nine[, c('dataset', 'age', 'gene', 'tested', 'logFC', 'CI.L', 'CI.R', 'P.Value')])
cat(readLines(file.path(out, 'SHA256SUMS')), sep = '\n')

#%% debug
