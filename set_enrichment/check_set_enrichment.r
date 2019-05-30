library(dplyr)
library(magrittr)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('f', 'config', 'yaml', '../config.yaml'),
    opt('c', 'cohort', 'chr', 'ACC'),
    opt('s', 'sets', 'RData', '../data/focal.RData'),
    opt('g', 'ng', 'RData', '../data/ng.RData'),
    opt('m', 'method', 'aracne or genenet', 'aracne'),
    opt('i', 'min_edges', 'integer', '1e3'),
    opt('a', 'max_edges', 'integer', '1e7'),
    opt('t', 'n_steps', 'integer', '15'),
    opt('n', 'network', 'RData', '../networks/aracne/naive/ACC.RData'),
    opt('o', 'outfile', '.RData', 'focal_aracne/naive/ACC.RData'))

config = io$read_yaml(args$config)
net = io$load(args$network)
ng = io$load(args$ng) %>% filter(cohort == args$cohort, expr == "naive") #TODO: also copycor
min_edges = as.numeric(args$min_edges)
max_edges = as.numeric(args$max_edges)
n_steps = as.numeric(args$n_steps)

valid_genes = valid_tfs = ng$genes[[1]]
if (args$method %in% config$has_tf)
    valid_tfs = ng$tfs[[1]]

set2possible_links = function(genes, net=NA) {
    ng = length(genes)
    ntf = sum(genes %in% valid_tfs)
    (ng - ntf) * ntf + 0.5 * (ntf - 1) * ntf
}
set2real_links = function(genes, net)
    nrow(filter(net, Regulator %in% genes & Target %in% genes))

# get gene sets of co-amplified segments (either focal CNA or aneuploidy)
rset = io$load(args$sets)
if (args$cohort %in% names(rset$sets)) { # separate sets for each cohort
    cna_genes = rset$sets[[args$cohort]]
} else { # common sets (e.g. chromosomes)
    cna_genes = rset$sets
}
cna_genes = gset$filter(cna_genes, min=2, max=Inf, valid=valid_genes)

# compare within-segment vs. outside of segment
#   and TF targets vs non-TF targets
calc_net = function(net, co) {
    subnet = head(net, co)
    do_test = function(seg_real, all_real, seg_psbl, all_psbl, ...) {
        links = rbind(c(seg_real, all_real), c(seg_psbl, all_psbl))
        colnames(links) = c("seg", "all")
        rownames(links) = c("real", "psbl")
        links[,2] = links[,2] - links[,1] # adjust all to outside of segment
        links[2,] = links[2,] - links[1,] # adjust possible to "not chosen"
        broom::tidy(fisher.test(links)) %>%
            select(-method, -alternative) %>%
            mutate(links = list(links))
    }
    make_counts = function(fun) {
        segs = sapply(cna_genes, function(g) fun(g, net=subnet))
        if (length(Reduce(intersect, cna_genes)) > 0)
            segs
        else # only merge if non-overlapping sets
            c(all=sum(segs), segs)
    }
    data.frame(seg_psbl = make_counts(set2possible_links), # within-segment links possible
               seg_real = make_counts(set2real_links), # within-segment links in network
               all_psbl = set2possible_links(valid_genes), # total links possible
               all_real = nrow(subnet)) %>% # total links in network
        tibble::rownames_to_column("seg_id") %>%
        as_tibble() %>%
        filter(seg_psbl > 0) %>%
        mutate(fet = purrr::pmap(., do_test)) %>%
        tidyr::unnest()
}
cutoff = exp(seq(log(min_edges), log(max_edges), length.out=n_steps))
if (nrow(net) < cutoff[length(cutoff)]) {
    cutoff = cutoff[cutoff <= nrow(net)]
    cutoff = exp(seq(log(min_edges), log(nrow(net)), length.out=length(cutoff)))
}
res = lapply(cutoff, calc_net, net=net) %>%
    setNames(cutoff) %>%
    bind_rows(.id="edges") %>%
    mutate(edges = as.integer(edges))

save(res, file=args$outfile)
