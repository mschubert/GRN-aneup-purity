sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')
aracne = import('tools/aracne')

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('m', 'method', 'method identifier', 'aracne'),
    opt('o', 'outfile', '.RData to save to', 'aracne/ACC.RData'))

expr = tcga$rna_seq(args$cohort, trans="vst")
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
expr = expr[!is.na(rownames(expr)) & rownames(expr) != "",]
net = aracne$aracne(expr)

save(net, file=args$outfile)
