library(dplyr)
library(magrittr)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('c', 'cohort', 'chr', 'ACC'),
    opt('r', 'regions', 'RData', '../data/focal.RData'),
    opt('m', 'method', 'aracne or genenet', 'aracne'),
    opt('n', 'network', 'RData', '../networks/aracne/ACC.RData'),
    opt('o', 'outfile', '.RData', 'focal_aracne/ACC.RData'))

util = import(paste0("./", args$method))
net = io$load(args$network)
valid_genes = util$valid_genes(net)

# get gene sets of co-amplified segments (either focal CNA or aneuploidy)
rset = io$load(args$regions)
keep = rownames(rset$copies[[args$cohort]])
cna_genes = rset$genes[keep] %>%
    gset$filter(min=2, max=Inf, valid=valid_genes)

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
