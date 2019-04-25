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

net = io$load(args$network)

if (args$method %in% c("aracne", "TFbinding")) { # with knowledge what is TF
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
#   and TF targets vs non-TF targets
calc_net = function(net, co) {
    net2 = head(net, co)
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
        segs = sapply(cna_genes, fun, net=net2)
        if (length(Reduce(intersect, cna_genes)) > 0)
            segs
        else # only merge if non-overlapping sets
            c(all=sum(segs), segs)
    }
    res = data.frame(seg_psbl = make_counts(set2possible_links), # filtered network
                     seg_real = make_counts(set2real_links),
                     all_psbl = set2possible_links(valid_genes, net), # full network
                     all_real = nrow(net)) %>%
        tibble::rownames_to_column("seg_id") %>%
        as_tibble() %>%
        filter(seg_psbl > 0) %>%
        mutate(fet = purrr::pmap(., do_test)) %>%
        tidyr::unnest()
}
cutoff = exp(seq(log(1e3), log(1e6), length.out=15))
if (nrow(net) < cutoff[length(cutoff)]) {
    cutoff = cutoff[cutoff <= nrow(net)]
    cutoff = exp(seq(log(1e3), log(nrow(net)), length.out=length(cutoff)))
}
res = lapply(cutoff, calc_net, net=net) %>%
    setNames(cutoff) %>%
    bind_rows(.id="edges") %>%
    mutate(edges = as.integer(edges))

save(res, file=args$outfile)
