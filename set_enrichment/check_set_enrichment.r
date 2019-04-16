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
    opt('n', 'network', 'RData', '../networks/aracne/ACC.RData'),
    opt('o', 'outfile', '.RData', 'focal_aracne/ACC.RData'))

util = import(paste0("./", args$method))
net = io$load(args$network)
valid_genes = util$valid_genes(net)

if (method == "aracne") {
    valid_genes = function(net) unique(c(net$Regulator, net$Target))
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
} else {
    valid_genes = function(net) unique(c(as.character(net$node1),
                                         as.character(net$node2)))
    set2possible_links = function(genes, net) {
        ng = length(genes)
        0.5 * (ng^2 - ng)
    }
    set2real_links = function(genes, net) {
        net %>%
            filter(node1 %in% genes & node2 %in% genes) %>%
            nrow()
    }
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
psbl = c(seg = sum(sapply(cna_genes, util$set2possible_links, net=net)),
         all = util$set2possible_links(valid_genes, net))
real = c(seg = sum(sapply(cna_genes, util$set2real_links, net=net)),
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
