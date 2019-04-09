library(dplyr)
sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('m', 'method', 'method identifier', 'aracne'),
    opt('o', 'outfile', '.RData to save to', 'aracne/ACC.RData'))

counts = tcga$rna_seq(args$cohort, trans="raw")
keep = rowSums(counts) >= 5 * ncol(counts)
rm(counts)

expr = tcga$rna_seq(args$cohort, trans="vst")[keep,]
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
expr = expr[!is.na(rownames(expr)) & rownames(expr) != "",]

switch(args$method,
    "aracne" = {
        gset = import('data/genesets')
        tfs = gset$go() %>%
            filter(id == "GO:0003700") %>%
            pull(hgnc_symbol)

        ar = import('tools/aracne')
        bs = 100
        clustermq::register_dopar_cmq(n_jobs=bs, memory=10240)
        net = ar$aracne(expr, tfs, folder=paste0(".temp_", args$cohort), bootstrap=bs)
    },
    "genenet" = {
        net = gnet$pcor(expr, fdr=0.05)
    },
    "TFbinding" = {
        sets = enr$genes("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
        names(sets) = sub("_[^_]+$", "", names(sets))

        net = stack(sets)[c(2,1)]
        colnames(net) = c("Regulator", "Target")
    },
    {
        fun = getFromNamespace(args$method, ns="netbenchmark")
        net = fun(t(expr)) %>%
            reshape2::melt() %>%
            filter(value != 0) %>%
            dplyr::rename(node1 = Var1, node2 = Var2)
    }
)

save(net, file=args$outfile)
