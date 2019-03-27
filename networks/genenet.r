sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')
gnet = import('tools/genenet')

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('m', 'method', 'method identifier', 'genenet'),
    opt('o', 'outfile', '.RData to save to', 'genenet/ACC.RData'))

counts = tcga$rna_seq(args$cohort, trans="raw")
keep = rowSums(counts) >= 5 * ncol(counts)
rm(counts)

expr = tcga$rna_seq(args$cohort, trans="vst")[keep,]
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
expr = expr[!is.na(rownames(expr)) & rownames(expr) != "",]
net = gnet$pcor(expr, fdr=0.05)

save(net, file=args$outfile)
