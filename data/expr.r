library(dplyr)
sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('r', 'correction', 'naive or copy-corrected', 'copycor'),
    opt('o', 'outfile', '.RData to save to', 'expr_copycor/ACC.RData'))

copies = na.omit(tcga$cna_genes(args$cohort, chr_excl=c("X", "Y", "MT")))
counts = tcga$rna_seq(args$cohort, trans="raw") %>%
    tcga$filter(cancer=TRUE, primary=TRUE)
counts = counts[rowSums(counts >= 10) > 0.2*ncol(counts),]
narray::intersect(counts, copies, along=1)
narray::intersect(counts, copies, along=2)

if (args$correction != "copycor")
    copies[] = 1

expr = DESeq2::DESeqDataSetFromMatrix(counts,
        colData = data.frame(id=colnames(counts)),
        design = ~1) %>%
    DESeq2::estimateSizeFactors(normMatrix=copies) %>%
    DESeq2::estimateDispersions(fitType="parametric") %>%
    DESeq2::getVarianceStabilizedData()

rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
expr = expr[!is.na(rownames(expr)) & !duplicated(rownames(expr)),]

save(expr, file=args$outfile)
