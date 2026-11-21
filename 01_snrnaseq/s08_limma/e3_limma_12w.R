# E3: S08 limma-voom with 12 weeks as a fifth age (Purkinje + Bergmann only).
# 2026-10-24. Declared extension; the frozen four-age release is not modified.
# Design: pseudobulk per mouse, TMM on the full retained transcriptome before filtering,
# limma-voom ~0+genotype_age, KI-WT at each of 5 ages, BH within the 1,414-gene panel per
# population x age; whole-transcriptome secondary BH within population x age.
# Usage: Rscript e3_limma_12w.R <FROZEN_PB_DIR> <NEW12_DIR> <CORE_REGISTRY> <FROZEN_RELEASE_HITS> <OUT_DIR>
# The inputs are external (see the README, "Data"); their as-run sha256 values are pinned below.
# Refuses to run if OUT_DIR exists or any pinned input hash differs.
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5L)
frozen <- args[1]; new12 <- args[2]; core_path <- args[3]; frozen_hits <- args[4]; out <- args[5]
if (dir.exists(out)) stop("output exists; refusing to overwrite: ", out)

PIN <- c(
  raw_counts = "3597cca1e4f690125b414ad48fe517c0f913d2f119f2741d118b43c451942778",
  unit_metadata = "11f96a5eaac14dcd523f65935bf255aa18caf2a734597dafee55919b3e97abb4",
  gene_metadata = "765cf0b8c4ae14711a87c3f64612b6144ca00451bfe6bc60d3fab7652ffa62a3",
  core_registry = "ea40dde47e744d9930df038458f4683264736d57e442d450cbc5f1d9707056b8",
  frozen_hits = "40593a55c8538dc7d687cc07eea13b8d640e4fe4d696794879d159e74b625298"
)
sha <- function(path) strsplit(system2("sha256sum", shQuote(path), stdout = TRUE), " ")[[1]][1]
stopifnot(identical(sha(file.path(frozen, "raw_counts.tsv.gz")), PIN[["raw_counts"]]))
stopifnot(identical(sha(file.path(frozen, "unit_metadata.tsv")), PIN[["unit_metadata"]]))
stopifnot(identical(sha(file.path(frozen, "gene_metadata.tsv")), PIN[["gene_metadata"]]))
stopifnot(identical(sha(core_path), PIN[["core_registry"]]))
stopifnot(identical(sha(frozen_hits), PIN[["frozen_hits"]]))

read_counts <- function(file) {
  r <- read.delim(gzfile(file), check.names = FALSE)
  ids <- r$stable_id; r$stable_id <- NULL
  m <- as.matrix(r); rownames(m) <- ids; m
}
prep <- function(counts) edgeR::calcNormFactors(edgeR::DGEList(counts = counts), method = "TMM")
expression_filter <- function(y, age) {
  stopifnot(length(age) == ncol(y))
  effective <- y$samples$lib.size * y$samples$norm.factors
  cutoff <- 10 / median(effective) * 1e6
  above <- edgeR::cpm(y, normalized.lib.sizes = TRUE) >= cutoff
  eligible <- rep(FALSE, nrow(y))
  for (a in sort(unique(age))) eligible <- eligible | rowSums(above[, age == a, drop = FALSE]) >= 3
  eligible & rowSums(y$counts) >= 15
}

fmeta <- read.delim(file.path(frozen, "unit_metadata.tsv"), colClasses = c(animal_id = "character"))
fmeta$genotype <- ifelse(fmeta$genotype == "SCA1", "KI", fmeta$genotype)
# Hard stop: the 24 discovery mice must equal the corrected labels (s06 file sha d2a0b8b2...).
corrected <- data.frame(
  animal_id = c("72931","72932","7294","7202","72921","72922",
                "6571","65731","65732","6569","65701","65702",
                "1589","2021","20212","1974","2020","20202",
                "58231","58241","58242","5812","5822","58232"),
  age_weeks = c(rep(5L,6), rep(18L,6), rep(24L,6), rep(30L,6)),
  genotype = c(rep(c("KI","KI","KI","WT","WT","WT"), 4)),
  sex = c("M","M","F","F","M","M", rep("M",6), "M","F","F","M","F","F", rep("F",6)),
  stringsAsFactors = FALSE
)
got_lab <- unique(fmeta[, c("animal_id", "age_weeks", "genotype", "sex")])
got_lab <- got_lab[order(got_lab$animal_id), ]
exp_lab <- corrected[order(corrected$animal_id), ]
stopifnot(nrow(got_lab) == 24L,
          identical(got_lab$animal_id, exp_lab$animal_id),
          identical(as.integer(got_lab$age_weeks), as.integer(exp_lab$age_weeks)),
          identical(as.character(got_lab$genotype), exp_lab$genotype),
          identical(as.character(got_lab$sex), exp_lab$sex))
genes_meta <- read.delim(file.path(frozen, "gene_metadata.tsv"), stringsAsFactors = FALSE)
core <- read.delim(core_path, stringsAsFactors = FALSE)
core <- core[tolower(as.character(core$primary_test_eligible_before_expression_filter)) == "true", ]
stopifnot(nrow(core) == 1414L)

fh <- read.delim(gzfile(frozen_hits), stringsAsFactors = FALSE)
fh_key <- ifelse(fh$population == "purkinje", "Purkinje cell",
          ifelse(fh$population == "bergmann", "Bergmann glia", NA_character_))
stopifnot(!anyNA(fh_key))
fh_hits <- aggregate(fh$decision == "HIT",
  list(population = fh_key, age_weeks = fh$age_weeks), sum)
names(fh_hits)[3] <- "frozen_core_hits"

new_meta <- read.delim(file.path(new12, "12w_unit_metadata.tsv"), colClasses = c(animal_id = "character"))
pops <- list(`Purkinje cell` = "12w_purkinje_raw_counts.tsv.gz",
             `Bergmann glia` = "12w_bergmann_raw_counts.tsv.gz")
tested_genes <- c("Kcnc3", "Kcnma1", "Kcnd3", "Cacna1a", "Cacna1g", "Itpr1",
                  "Atp2a3", "Scn8a", "Hcn1", "Kcnn2", "Kcna1")
ages_expected <- c(5L, 12L, 18L, 24L, 30L)

all_rows <- list(); hits <- list(); gene_rows <- list(); checks <- list()
for (pop in names(pops)) {
  fc <- read_counts(file.path(frozen, "raw_counts.tsv.gz"))
  fm <- fmeta[fmeta$cell_type_raw == pop, ]
  stopifnot(nrow(fm) == 24L, all(fm$unit_label %in% colnames(fc)))
  fc <- fc[, fm$unit_label, drop = FALSE]
  nc <- read_counts(file.path(new12, pops[[pop]]))
  nm <- new_meta[new_meta$cell_type == pop, ]
  nm <- nm[match(colnames(nc), nm$unit_label), ]
  stopifnot(nrow(nm) == 6L, all(nm$age_weeks == 12L))
  if (!setequal(rownames(fc), rownames(nc))) stop("gene ID sets differ between frozen and 12w counts")
  nc <- nc[rownames(fc), , drop = FALSE]
  counts <- cbind(fc, nc)
  meta <- rbind(
    data.frame(animal_id = fm$animal_id, age_weeks = as.integer(fm$age_weeks),
               genotype = fm$genotype, sex = fm$sex, unit_label = fm$unit_label,
               stringsAsFactors = FALSE),
    data.frame(animal_id = nm$animal_id, age_weeks = 12L, genotype = nm$genotype,
               sex = nm$sex, unit_label = nm$unit_label, stringsAsFactors = FALSE))
  stopifnot(ncol(counts) == 30L, !anyDuplicated(meta$animal_id))
  stopifnot(identical(as.integer(sort(unique(meta$age_weeks))), as.integer(ages_expected)))
  stopifnot(all(table(meta$age_weeks, meta$genotype) == 3L))

  y <- prep(counts)
  full_sizes <- y$samples$lib.size; nf <- y$samples$norm.factors
  keep <- expression_filter(y, meta$age_weeks)
  y <- y[keep, , keep.lib.sizes = TRUE]
  stopifnot(all(y$samples$lib.size == full_sizes), all(y$samples$norm.factors == nf))
  core_ids <- core$stable_id[core$stable_id %in% rownames(y)]
  core_idx <- match(core_ids, rownames(y))
  stopifnot(length(core_idx) > 0L)  # expression filter removes panel genes; tested set is recorded

  genotype <- factor(meta$genotype, levels = c("WT", "KI"))
  agef <- factor(meta$age_weeks, levels = ages_expected)
  genotype_age <- interaction(genotype, agef, sep = "_", lex.order = TRUE)
  design <- model.matrix(~0 + genotype_age)
  colnames(design) <- levels(genotype_age)
  col_for <- function(g, a) {
    hit <- grep(paste0("^", g, "_", a, "$"), colnames(design), value = TRUE)
    stopifnot(length(hit) == 1L); hit
  }
  contrasts <- matrix(0, ncol(design), length(ages_expected),
                      dimnames = list(colnames(design), paste0("age", ages_expected)))
  for (k in seq_along(ages_expected)) {
    contrasts[col_for("KI", ages_expected[k]), k] <- 1
    contrasts[col_for("WT", ages_expected[k]), k] <- -1
  }
  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v$E, design, weights = v$weights)
  fit <- eBayes(contrasts.fit(fit, contrasts), robust = TRUE)
  for (k in seq_along(ages_expected)) {
    tt <- topTable(fit, coef = k, number = Inf, confint = 0.95, sort.by = "none")
    tt$stable_id <- rownames(tt)
    tt$population <- pop; tt$age_weeks <- ages_expected[k]
    tt$dataset_symbol <- genes_meta$dataset_symbol[match(tt$stable_id, genes_meta$stable_id)]
    finite <- is.finite(tt$logFC) & is.finite(tt$P.Value) & is.finite(tt$CI.L) & is.finite(tt$CI.R)
    p <- ifelse(finite, tt$P.Value, 1)
    core_hit <- tt$stable_id %in% core_ids
    q_core <- rep(NA_real_, nrow(tt)); q_core[core_hit] <- p.adjust(p[core_hit], "BH")
    q_sec <- p.adjust(p, "BH")
    frac_unresolved <- mean(!finite[core_hit])
    if (frac_unresolved > 0.01)
      stop(sprintf("STOP RULE: >1%% unresolved core tests in %s age %d (%.4f)",
                   pop, ages_expected[k], frac_unresolved))
    checks[[length(checks) + 1L]] <- data.frame(population = pop, age_weeks = ages_expected[k],
      core_tests = sum(core_hit), core_unresolved = sum(!finite[core_hit]),
      core_BH_hits = sum(q_core < 0.05, na.rm = TRUE),
      frozen_four_age_core_hits = fh_hits$frozen_core_hits[
        match(paste(pop, ages_expected[k]), paste(fh_hits$population, fh_hits$age_weeks))],
      secondary_tests = nrow(tt), secondary_BH_hits = sum(q_sec < 0.05),
      fitted_genes = nrow(y), filtered_out = nrow(counts) - nrow(y))
    all_rows[[length(all_rows) + 1L]] <- data.frame(
      stable_id = tt$stable_id, dataset_symbol = tt$dataset_symbol, population = pop,
      age_weeks = ages_expected[k], log2FC = tt$logFC, CI_L = tt$CI.L, CI_R = tt$CI.R,
      p = p, q_core = q_core, q_secondary = q_sec,
      core_panel = core_hit, failed_fit_p1 = !finite,
      status_core = ifelse(q_core < 0.05, "HIT", "TESTED, NOT A HIT"),
      status_secondary = ifelse(q_sec < 0.05, "HIT", "TESTED, NOT A HIT"),
      stringsAsFactors = FALSE)
    if (ages_expected[k] == 12L) {
      g <- all_rows[[length(all_rows)]]
      gene_rows[[pop]] <- g[match(tested_genes, g$dataset_symbol),
                            c("dataset_symbol", "population", "age_weeks", "log2FC",
                              "CI_L", "CI_R", "p", "q_core", "q_secondary", "status_core")]
    }
  }
  cat("population", pop, "fitted genes", nrow(y), "core present", length(core_idx), "\n")
}
res <- do.call(rbind, all_rows)
dir.create(out, recursive = TRUE)
write.table(res, file.path(out, "LIMMA_RESULTS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
chk <- do.call(rbind, checks)
write.table(chk, file.path(out, "HITS_BY_AGE.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
tw <- do.call(rbind, gene_rows)
rownames(tw) <- NULL
write.table(tw, file.path(out, "TWELVE_WEEK_GENES.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))

L <- c("# S08 limma-voom with 12 weeks as a fifth age (declared extension)", "",
  "Purkinje cells and Bergmann glia only. The frozen four-age release is unchanged.", "",
  "## Hits (BH < 0.05)", "",
  "| population | age | core panel hits | frozen four-age core hits | secondary hits | core tests | fitted genes |",
  "|---|---:|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(chk))) L <- c(L, sprintf("| %s | %d | %d | %s | %d | %d | %d |",
  chk$population[i], chk$age_weeks[i], chk$core_BH_hits[i],
  ifelse(is.na(chk$frozen_four_age_core_hits[i]), "n/a", as.character(chk$frozen_four_age_core_hits[i])),
  chk$secondary_BH_hits[i], chk$core_tests[i], chk$fitted_genes[i]))
L <- c(L, "", "## 12-week values for the requested genes", "",
  "| gene | population | log2FC | 95% CI | p | q (1,414 panel) | q (secondary) |", "|---|---|---:|---|---:|---:|---:|")
for (i in seq_len(nrow(tw))) L <- c(L, sprintf("| %s | %s | %.3f | %.3f to %.3f | %.3g | %.3g | %.3g |",
  tw$dataset_symbol[i], tw$population[i], tw$log2FC[i], tw$CI_L[i], tw$CI_R[i],
  tw$p[i], tw$q_core[i], tw$q_secondary[i]))
writeLines(L, file.path(out, "SUMMARY.md"))
rec <- list(stage = "E3 S08 12-week extension", created_utc = format(Sys.time(), tz = "UTC"),
  script_sha256 = sha(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
  frozen = frozen, new12 = new12, core_registry = core_path,
  frozen_release_hits = frozen_hits, frozen_release_hits_sha256 = PIN[["frozen_hits"]],
  discovery_labels = "24 discovery mice match the corrected s06 labels (unit_metadata sha d2a0b8b2...)",
  sex_note = paste("The six 12-week mice are all male (animal manifest sha 7c2f2d4d..., checked in E5);",
                   "there is no within-age sex contrast at 12 weeks, so no sex term is used."),
  h5ad_note = "12w counts built from the pinned annotated h5ad by the 12-week pseudobulk preparation script (not in this package)",
  ages = ages_expected, populations = names(pops), design = "~0+genotype_age; KI-WT per age",
  normalization = "TMM on full retained transcriptome before filtering; keep.lib.sizes TRUE",
  family = "BH within the 1,414-gene panel per population x age; secondary BH within population x age",
  stop_rule = "core unresolved > 1% stops the run", hits = chk,
  outputs = setNames(vapply(c("LIMMA_RESULTS.tsv", "HITS_BY_AGE.tsv", "TWELVE_WEEK_GENES.tsv", "SUMMARY.md"),
    function(n) sha(file.path(out, n)), ""), c("LIMMA_RESULTS.tsv", "HITS_BY_AGE.tsv", "TWELVE_WEEK_GENES.tsv", "SUMMARY.md")))
write_json(rec, file.path(out, "RECORD.json"), auto_unbox = TRUE, pretty = TRUE)
cat(readLines(file.path(out, "SUMMARY.md")), sep = "\n")

#%% debug
