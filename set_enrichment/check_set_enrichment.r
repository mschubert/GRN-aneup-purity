library(dplyr)
library(magrittr)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('c', 'cohort', 'chr', 'ACC'),
    opt('s', 'sets', 'RData', '../data/focal.RData'),
    opt('m', 'method', 'aracne or genenet', 'aracne'),
    opt('n', 'network', 'RData', '../networks/aracne/naive/ACC.RData'),
    opt('o', 'outfile', '.RData', 'focal_aracne/naive/ACC.RData'))

util = import(paste0("./", args$method))
net = io$load(args$network)

if (args$method == "aracne") { # for now, only one with knowledge what is TF
    valid_genes = unique(c(net$Regulator, net$Target))
    set2possible_links = function(genes, net) {
        ntf = sum(genes %in% net$Regulator)
        ntg = sum(genes %in% setdiff(net$Target, net$Regulator))
        if (ntf > 0)
            ntf * (ntf-1) + ntf * ntg
        else
            0
    }
    set2real_links = function(genes, net)
        nrow(filter(net, Regulator %in% genes & Target %in% genes))
} else {
    valid_genes = unique(c(as.character(net$node1), as.character(net$node2)))
    set2possible_links = function(genes, net) { ng = length(genes); 0.5 * (ng^2 - ng) }
    set2real_links = function(genes, net)
        nrow(filter(net, node1 %in% genes & node2 %in% genes))
}

# get gene sets of co-amplified segments (either focal CNA or aneuploidy)
rset = io$load(args$sets)
if (args$cohort %in% names(rset$sets)) { # separate sets for each cohort
    cna_genes = rset$sets[[args$cohort]]
} else { # common sets (e.g. chromosomes)
    cna_genes = rset$sets
}
cna_genes = gset$filter(cna_genes, min=2, max=Inf, valid=valid_genes)

# compare within-segment vs. outside of segment
#   and TF targets vs non-TF targets (FET?)
#TODO: also check for each segment individually
psbl = c(seg = sum(sapply(cna_genes, set2possible_links, net=net)),
         all = set2possible_links(valid_genes, net))
real = c(seg = sum(sapply(cna_genes, set2real_links, net=net)),
         all = nrow(net))
links = rbind(real, psbl) # so estimate is OR segments over rest
links[,2] = links[,2] - links[,1] # adjust all to outside of segment
links[2,] = links[2,] - links[1,] # adjust possible to "not chosen"
res = broom::tidy(fisher.test(links)) %>%
    select(-method) %>%
    mutate(links = list(links))

res
save(res, file=args$outfile)

# plot total number of amp segments (/genes?) per cancer type
# plot enrichment of amp/del in segments per method
