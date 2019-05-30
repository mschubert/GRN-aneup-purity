library(dplyr)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')

#' Convert a network matrix to a data.frame listing edges
#'
#' @param mat  link scores, high is high confidence (i.e., not p-values)
#' @return  A data.frame with fields: Regulator, Target, score
mat2df = function(mat, tfs) {
    # symmetrise regulators only (or all genes if no tfs)
    mat[tfs, tfs] = pmax(mat[tfs, tfs], t(mat[tfs, tfs]))
    # make sure we don't have the same links twice
    lvls = unique(c(rownames(mat), colnames(mat)))
    reshape2::melt(mat) %>%
        transmute(Regulator = factor(Var1, levels=lvls),
                  Target = factor(Var2, levels=lvls),
                  score = value) %>%
        filter(score != 0,
               as.integer(Regulator) < as.integer(Target)) %>%
        arrange(-score) %>%
        mutate(Regulator = as.character(Regulator),
               Target = as.character(Target))
}

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('e', 'expr', 'expression matrix RData', '../data/expr_copycor/ACC.RData'),
    opt('s', 'select', 'top N links max', '1e7'),
    opt('m', 'method', 'method identifier', 'aracne'),
    opt('a', 'tf_annot', 'RData', '../data/tf_annot.RData'),
    opt('b', 'tf_binding', 'RData', '../data/tf_binding.RData'),
    opt('o', 'outfile', '.RData to save to', 'aracne/copycor/ACC.RData'))

expr = io$load(args$expr)
gs = rownames(expr)
tfs = intersect(io$load(args$tf_annot), rownames(expr))
top_n = as.numeric(args$select)

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
            arrange(qval, pval) %>%
            transmute(Regulator = as.character(node1),
                      Target = as.character(node2),
                      pval = pval, qval=qval)
    },
    "Genie3" = {
        clustermq::register_dopar_cmq(n_jobs=200, memory=2048)
        net = GENIE3::GENIE3(expr, verbose=TRUE) %>%
            mat2df(tfs=gs)
    },
    "Genie3+TF" = {
        clustermq::register_dopar_cmq(n_jobs=100, memory=2048)
        net = GENIE3::GENIE3(expr, regulators=tfs, verbose=TRUE) %>%
            mat2df(tfs=tfs)
    },
    "Tigress" = {
        clustermq::register_dopar_cmq(n_jobs=500, memory=8192)
        net = tigress::tigress(t(expr), allsteps=FALSE, verb=TRUE,
                   usemulticore=TRUE) %>% mat2df(tfs=gs)
    },
    "Tigress+TF" = {
        clustermq::register_dopar_cmq(n_jobs=200, memory=2048)
        net = tigress::tigress(t(expr), tflist=tfs, allsteps=FALSE,
                   verb=TRUE, usemulticore=TRUE) %>% mat2df(tfs=tfs)
    },
    "TFbinding" = {
        sets = io$load(args$tf_binding)
        net = stack(sets)[c(2,1)]
        colnames(net) = c("Regulator", "Target")
    },
    {
        fun = getFromNamespace(args$method, ns="netbenchmark")
        net = fun(t(expr)) %>% mat2df(tfs=gs)
    }
)

net = head(net, top_n)
save(net, file=args$outfile)
