sys = import('sys')
enr = import('tools/enrichr')

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', ''),
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('o', 'outfile', 'file to save to', 'tfs.RData'))

sets = enr$genes("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
names(sets) = sub("_[^_]+$", "", names(sets))

save(sets, file=args$outfile)
