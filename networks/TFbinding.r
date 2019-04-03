sys = import('sys')
enr = import('tools/enrichr')

args = sys$cmd$parse(
    opt('c', 'cohort', 'TCGA identifier (ignored)', 'ACC'),
    opt('m', 'method', 'method name (ignored)', 'NA'),
    opt('b', 'bootstraps', 'network inference (ignored)', 'NA'),
    opt('o', 'outfile', 'file to save to', 'TFbinding/ACC.RData'))

sets = enr$genes("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
names(sets) = sub("_[^_]+$", "", names(sets))

net = stack(sets)[c(2,1)]
colnames(net) = c("Regulator", "Target")

save(net, file=args$outfile)
