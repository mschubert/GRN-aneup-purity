library(dplyr)
sys = import('sys')
gset = import('data/genesets')
enr = import('tools/enrichr')

args = sys$cmd$parse(
    opt('a', 'annot', 'output RData', 'tf_annot.RData'),
    opt('b', 'binding', 'output RData', 'tf_binding.RData'))

tf_annot = gset$go() %>%
    filter(id == "GO:0003700") %>%
    pull(hgnc_symbol)
save(tf_annot, file=args$annot)

tf_binding = enr$genes("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
names(tf_binding) = sub("_[^_]+$", "", names(tf_binding))
save(tf_binding, file=args$binding)
