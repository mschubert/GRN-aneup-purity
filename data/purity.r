io = import('io')
sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')

cohort2genes = function(cohort) {
    idx = as.data.frame(estimate[[cohort]])
    eset = tcga$rna_seq(cohort) %>%
        tcga$filter(primary=TRUE, cancer=TRUE)
    common = intersect(rownames(idx), colnames(eset))
    eset = eset[,common] %>%
        DESeq2::DESeqDataSetFromMatrix(colData=idx[common,,drop=FALSE], ~1+purity) %>%
        DESeq2::estimateSizeFactors()
    rownames(eset) = idmap$gene(rownames(eset), to="external_gene_name")
    eset = eset[rowMeans(counts(eset)) >= 5 &
                !is.na(rownames(eset)) &
                !duplicated(rownames(eset)),]
    diff_expr = DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomLRT(reduced=~1, maxit=1000) %>%
        DESeq2::results() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("hgnc_symbol") %>%
        dplyr::arrange(padj)
    cuts = c(100, 500, 1000, 2000, 3000, 4000, 5000)
    genes = lapply(cuts, function(n) head(diff_expr$hgnc_symbol, n)) %>%
        setNames(paste0("purity", cuts))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', ''),
    opt('o', 'outfile', 'file to save to', 'purity.RData'))

purity = tcga$purity() %>%
    filter(!is.na(estimate))

estimate = matrix(purity$estimate, ncol=1, dimnames=list(purity$Sample, "purity")) %>%
    narray::split(along=1, subsets=purity$cohort)
sets = clustermq::Q(cohort2genes, names(estimate), job_size=1, memory=8192,
                     export=list(tcga=tcga, idmap=idmap, estimate=estimate)) %>%
    setNames(names(estimate))

save(estimate, sets, file=args$outfile)
