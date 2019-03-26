sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')
gset = import('data/genesets')
ar = import('tools/aracne')

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('m', 'method', 'method identifier', 'aracne'),
    opt('b', 'bootstraps', 'number of bootstraps', '100'),
    opt('o', 'outfile', '.RData to save to', 'aracne/ACC.RData'))

tfs = gset$go() %>%
    filter(id == "GO:0003700") %>%
    pull(hgnc_symbol)

expr = tcga$rna_seq(args$cohort, trans="vst")
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
expr = expr[!is.na(rownames(expr)) & rownames(expr) != "",]

bs = as.integer(args$bootstraps)
clustermq::register_dopar_cmq(n_jobs=bs, memory=10240)
net = ar$aracne(expr, tfs, folder=paste0(".temp_", args$cohort), bootstrap=bs)

save(net, file=args$outfile)
