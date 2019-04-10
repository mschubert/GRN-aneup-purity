library(dplyr)
sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('s', 'select', 'top N links max', '1e6'),
    opt('m', 'method', 'method identifier', 'aracne'),
    opt('o', 'outfile', '.RData to save to', 'aracne/ACC.RData'))

counts = tcga$rna_seq(args$cohort, trans="raw")
keep = rowSums(counts) >= 5 * ncol(counts)
rm(counts)

expr = tcga$rna_seq(args$cohort, trans="vst")[keep,]
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
expr = expr[!is.na(rownames(expr)) & rownames(expr) != "",]

top_n = as.integer(args$select)

switch(args$method,
    "aracne" = {
        gset = import('data/genesets')
        tfs = gset$go() %>%
            filter(id == "GO:0003700") %>%
            pull(hgnc_symbol)

        ar = import('tools/aracne')
        bs = 100
        clustermq::register_dopar_cmq(n_jobs=bs, memory=10240)
        net = ar$aracne(expr, tfs, folder=paste0(".temp_", args$cohort),
                        bootstrap=bs, p.value=0.1) %>%
            arrange(pvalue) %>%
            head(top_n)
    },
    "genenet" = {
        net = gnet$pcor(expr, fdr=1) %>%
            arrange(qval, pval) %>%
            head(top_n)
    },
    "TFbinding" = {
        sets = enr$genes("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
        names(sets) = sub("_[^_]+$", "", names(sets))

        net = stack(sets)[c(2,1)]
        colnames(net) = c("Regulator", "Target")
    },
    {
        fun = getFromNamespace(args$method, ns="netbenchmark")
        net = fun(t(expr))
        net = reshape2::melt(net) %>%
            mutate(node1 = factor(Var1, levels=colnames(net)),
                   node2 = factor(Var2, levels=colnames(net)),
                   score = value) %>%
            filter(value != 0, node1 < node2) %>%
            arrange(-value) %>%
            head(top_n)
    }
)

save(net, file=args$outfile)
