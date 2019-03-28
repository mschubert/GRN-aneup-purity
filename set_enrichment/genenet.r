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
