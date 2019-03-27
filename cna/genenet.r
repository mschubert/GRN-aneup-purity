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
