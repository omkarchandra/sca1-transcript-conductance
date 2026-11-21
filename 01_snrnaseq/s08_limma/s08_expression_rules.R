# Label-invariant filtering and full-transcriptome normalization for frozen S01.
prepare_pseudobulk <- function(counts) {
    y <- edgeR::DGEList(counts=counts)
    edgeR::calcNormFactors(y,method="TMM")
}

expression_filter <- function(y, age) {
    stopifnot(length(age)==ncol(y))
    effective <- y$samples$lib.size*y$samples$norm.factors
    cutoff <- 10/median(effective)*1e6
    above <- edgeR::cpm(y,normalized.lib.sizes=TRUE)>=cutoff
    eligible <- rep(FALSE,nrow(y))
    for (a in sort(unique(age))) eligible <- eligible | rowSums(above[,age==a,drop=FALSE])>=3
    eligible & rowSums(y$counts)>=15
}

#%% debug
