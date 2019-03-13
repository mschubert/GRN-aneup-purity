sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')
gnet = import('tools/genenet')

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('m', 'method', 'method identifier', 'genenet'),
    opt('o', 'outfile', '.RData to save to', 'genenet/ACC.RData'))

expr = tcga$rna_seq(args$cohort, trans="vst")
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
expr = expr[!is.na(rownames(expr)) & rownames(expr) != "",]
net = gnet$pcor(expr)

save(net, file=args$outfile)
