library(dplyr)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('e', 'expr', 'expression matrix RData', '../data/expr_copycor/ACC.RData'),
    opt('s', 'select', 'top N links max', '1e6'),
    opt('m', 'method', 'method identifier', 'aracne'),
    opt('a', 'tf_annot', 'RData', '../data/tf_annot.RData'),
    opt('b', 'tf_binding', 'RData', '../data/tf_binding.RData'),
    opt('o', 'outfile', '.RData to save to', 'aracne/copycor/ACC.RData'),
    opt('p', 'parallel', 'number of cores', '1'))

expr = io$load(args$expr)
tfs = io$load(args$tf_annot)
cores = as.integer(args$parallel)

switch(args$method,
    "aracne" = {
        ar = import('tools/aracne')
        bs = 100
        clustermq::register_dopar_cmq(n_jobs=bs, memory=10240)
        temp = paste0(".temp", args$cohort,
                      paste(sample(letters, 5), collapse=""), sep="_")
        net = ar$aracne(expr, tfs, folder=temp, bootstrap=bs) %>%
            arrange(pvalue)
    },
    "genenet" = {
        gnet = import('tools/genenet')
        net = gnet$pcor(expr, fdr=0.01) %>%
            arrange(qval, pval)
    },
    "Genie3" = {
        net = GENIE3::GENIE3(expr, nCores=cores) %>%
            GENIE3::getLinkList() %>%
            arrange(-weight) %>%
            transmute(Regulator = regulatoryGene,
                      Target = targetGene,
                      score = weight)
    },
    "Genie3+TF" = {
        net = GENIE3::GENIE3(expr, regulators=tfs, nCores=cores) %>%
            GENIE3::getLinkList() %>%
            arrange(-weight) %>%
            transmute(Regulator = regulatoryGene,
                      Target = targetGene,
                      score = weight)
    },
    "Tigress" = {
        net = tigress::tigress(t(expr), usemulticore=cores)
    },
    "Tigress+TF" = {
        net = tigress::tigress(t(expr), tflist=tfs, usemulticore=cores)
    },
    "TFbinding" = {
        sets = io$load(args$tf_binding)
        net = stack(sets)[c(2,1)]
        colnames(net) = c("Regulator", "Target")
    },
    {
        fun = getFromNamespace(args$method, ns="netbenchmark")
        net = fun(t(expr))
        net = reshape2::melt(net) %>%
            transmute(node1 = factor(Var1, levels=colnames(net)),
                      node2 = factor(Var2, levels=colnames(net)),
                      score = value) %>%
            filter(score != 0, as.integer(node1) < as.integer(node2)) %>%
            arrange(-score)
    }
)

net = head(net, as.integer(args$select))
save(net, file=args$outfile)
