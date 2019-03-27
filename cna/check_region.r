library(dplyr)
library(magrittr)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('m', 'meta', 'RData', '../data/meta.RData'),
    opt('t', 'tfs', 'RData', '../data/tfs.RData'),
    opt('c', 'cohort', 'chr', 'ACC'),
#    opt('r', 'ranges', 'RData', '../'),
    opt('n', 'network', 'RData', '../networks/aracne/ACC.RData'),
    opt('o', 'outfile', '.RData', 'ACC.RData'))

meta = io$load(args$meta)
tfs = io$load(args$tfs)
net = io$load(args$network)

focal = meta$focal %>%
    filter(tcga$barcode2study(Sample) == args$cohort)

net_all = unique(c(net$Regulator, net$Target))
net_tfs = unique(net$Regulator)
net_tgs = setdiff(net_all, net_tfs)

# quantify the copy number of every tissue-specific segment
#   divide this by purity to get the cancer copies
cnas = narray::construct(copies ~ Sample + name, data=focal)
cna_genes = ungroup(focal) %>%
    select(name, gene) %>%
    filter(!duplicated(name)) %$%
    setNames(gene, name) %>%
    gset$filter(min=2, max=Inf, valid=net_all)

# compare within-segment vs. outside of segment
#   and TF targets vs non-TF targets (FET?)
set2possible_links = function(genes) {
    ntf = sum(genes %in% net_tfs)
    ntg = sum(genes %in% net_tgs)
    if (ntf > 0)
        #ntf-1 + ntg
        ntf * (ntf-1) + ntf * ntg
    else
        0
}
set2real_links = function(genes, net) {
    net %>%
        filter(Regulator %in% genes & Target %in% genes) %>%
        nrow()
}
possible = c(within = sum(sapply(cna_genes, set2possible_links)),
             outside = set2possible_links(net_all) - within_segs)
real = c(within = sum(sapply(cna_genes, set2real_links, net=net)),
         outside = nrow(net))
rbind(real, possible=possible-real) %>%
    fisher.test() %>%
    broom::tidy()

# plot total number of amp segments (/genes?) per cancer type

# plot enrichment of amp/del in segments per method
