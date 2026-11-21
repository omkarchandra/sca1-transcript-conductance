source(file.path(Sys.getenv("RELEASE_ROOT", getwd()), "config", "paths.R"))
# A2 step (a): per-gene negative-binomial overdispersion for the NB
# count split, estimated once on the full Purkinje counts. edgeR estimateDisp on the S08 Purkinje pseudobulk (24 discovery
# mice): genes with any count, TMM, design ~0 + genotype_age (as S08), robust = TRUE. NB size b = 1 / tagwise dispersion.
# Read-only; writes a NEW real_v01/dispersion/ next to this script.
#   Rscript ${DATA_ROOT}/publication/A2_nb_countsplit/a2_dispersion.R
suppressPackageStartupMessages(library(edgeR))
root <- release_expand('${REPO_ROOT}')
here <- file.path(root, release_expand('${DATA_ROOT}/publication/A2_nb_countsplit'))
out <- file.path(here, 'real_v01', 'dispersion')
pins <- c(counts = release_expand('${DATA_ROOT}/results/s08_pseudobulk_v01/raw_counts.tsv.gz'), units = release_expand('${DATA_ROOT}/results/s08_pseudobulk_v01/unit_metadata.tsv'))
want <- c(counts = '3597cca1e4f690125b414ad48fe517c0f913d2f119f2741d118b43c451942778',
          units = '11f96a5eaac14dcd523f65935bf255aa18caf2a734597dafee55919b3e97abb4')
sha <- function(p) strsplit(system2('sha256sum', shQuote(p), stdout = TRUE), ' ')[[1]][1]
for (k in names(pins)) if (!identical(sha(file.path(root, pins[[k]])), want[[k]])) stop('input hash differs: ', k)
if (dir.exists(file.path(here, 'real_v01'))) stop('real_v01 already exists')

raw <- read.delim(gzfile(file.path(root, pins[['counts']])), check.names = FALSE, stringsAsFactors = FALSE)
units <- read.delim(file.path(root, pins[['units']]), stringsAsFactors = FALSE, colClasses = c(animal_id = 'character'))
pk <- units[units$cell_type_raw == 'Purkinje cell', ]
stopifnot(nrow(pk) == 24, all(pk$unit_label %in% names(raw)))
cnt <- as.matrix(raw[, pk$unit_label]); rownames(cnt) <- raw$stable_id
any_count <- rowSums(cnt) > 0
y <- DGEList(counts = cnt[any_count, ])
y <- calcNormFactors(y, method = 'TMM')
group <- factor(paste(pk$genotype, pk$age_weeks, sep = '_'))
design <- model.matrix(~0 + group)
y <- estimateDisp(y, design, robust = TRUE)
tab <- data.frame(stable_id = rownames(y), total_count = rowSums(y$counts), AveLogCPM = y$AveLogCPM,
                  tagwise_dispersion = y$tagwise.dispersion, trended_dispersion = y$trended.dispersion,
                  nb_size_b = 1 / y$tagwise.dispersion, stringsAsFactors = FALSE)
dir.create(file.path(here, 'real_v01')); dir.create(out)
con <- file(file.path(out, 'DISPERSION.tsv'), open = 'wx')
writeLines('# A2: edgeR estimateDisp (robust) on the S08 Purkinje pseudobulk, 24 mice, TMM, ~0+genotype_age; b = 1/tagwise dispersion', con)
write.table(tab, con, sep = '\t', quote = FALSE, row.names = FALSE); close(con)
rec <- list(script_sha256 = sha(file.path(here, 'a2_dispersion.R')), inputs = as.list(setNames(want, pins)),
            n_genes_with_counts = sum(any_count), n_genes_zero = sum(!any_count), common_dispersion = y$common.dispersion,
            tagwise_quantiles = as.list(quantile(y$tagwise.dispersion, c(0, .1, .5, .9, 1))),
            design = '~0 + genotype_age', edgeR = as.character(packageVersion('edgeR')), R = R.version.string)
con <- file(file.path(out, 'RECORD.json'), open = 'wx'); writeLines(jsonlite::toJSON(rec, auto_unbox = TRUE, pretty = TRUE, digits = NA), con); close(con)
con <- file(file.path(out, 'SHA256SUMS'), open = 'wx')
for (f in sort(list.files(out))) if (f != 'SHA256SUMS') writeLines(paste0(sha(file.path(out, f)), '  ', f), con)
close(con)
print(rec[c('n_genes_with_counts', 'common_dispersion', 'tagwise_quantiles')])
cat(readLines(file.path(out, 'SHA256SUMS')), sep = '\n')

#%% debug
