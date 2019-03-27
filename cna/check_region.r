library(dplyr)
library(magrittr)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('c', 'cohort', 'chr', 'ACC'),
    opt('r', 'regions', 'RData', '../data/focal.RData'),
    opt('t', 'tfs', 'RData', '../data/tfs.RData'),
    opt('n', 'network', 'RData', '../networks/aracne/ACC.RData'),
    opt('o', 'outfile', '.RData', 'aracne/ACC.RData'))

tfs = io$load(args$tfs)

net = io$load(args$network)
valid_genes = unique(c(net$Regulator, net$Target))

# get gene sets of co-amplified segments (either focal CNA or aneuploidy)
rset = io$load(args$regions)
keep = rownames(rset$copies[[args$cohort]])
cna_genes = rset$genes[keep] %>%
    gset$filter(min=2, max=Inf, valid=valid_genes)

# compare within-segment vs. outside of segment
#   and TF targets vs non-TF targets (FET?)
set2possible_links = function(genes, net) {
    ntf = sum(genes %in% net$Regulator)
    ntg = sum(genes %in% setdiff(net$Target, net$Regulator))
    if (ntf > 0)
        ntf * (ntf-1) + ntf * ntg
    else
        0
}
set2real_links = function(genes, net) {
    net %>%
        filter(Regulator %in% genes & Target %in% genes) %>%
        nrow()
}
#TODO: also check for each segment individually
psbl = c(seg = sum(sapply(cna_genes, set2possible_links, net=net)),
         all = set2possible_links(valid_genes, net))
real = c(seg = sum(sapply(cna_genes, set2real_links, net=net)),
         all = nrow(net))
links = rbind(real, psbl) # so estimate is OR segments over rest
links[,2] = links[,2] - links[,1] # adjust all to outside of segment
links[2,] = links[2,] - links[1,] # adjust possible to "not chosen"
res = broom::tidy(fisher.test(links))

res
save(res, file=args$outfile)

# plot total number of amp segments (/genes?) per cancer type
# plot enrichment of amp/del in segments per method
