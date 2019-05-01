library(dplyr)
library(magrittr)
library(cowplot)
io = import('io')
sys = import('sys')
st = import('stats')

check_hits = function(net, real) {
    re = c(paste(as.character(net[[1]]), as.character(net[[2]])),
           paste(as.character(net[[2]]), as.character(net[[1]]))) %in%
        paste(as.character(real[[1]]), as.character(real[[2]])) + 0
    roc = data.frame(TP = c(0,cumsum(re)), FP = c(0,cumsum(1-re)))
    roc2 = roc[re != 0 | c(1,re[-c(1,2)], 1) != 0,]
}

args = sys$cmd$parse(
    opt('r', 'real', 'RData', '../networks/TFbinding/naive/ACC.RData'),
    opt('g', 'ng', 'number of genes per expr', '../data/ng.RData'),
    opt('o', 'outfile', '.RData', 'TFbinding_enrichment.RData'),
    opt('p', 'plotfile', 'pdf', 'TFbinding_enrichment.pdf'),
    arg('net', 'networks', arity='*', rep('../networks/aracne/naive/ACC.RData', 2)))

real = unique(io$load(args$real))
ng = io$load(args$ng) %>%
    mutate(slope = nrow(real)/psbl,
           slope_tf = nrow(real)/psbl_tf) %>%
    select(cohort, expr, slope, slope_tf) %>%
    tidyr::gather("kind", "slope", -cohort, -expr)

hits = do.call(rbind, strsplit(tools::file_path_sans_ext(args$net), "/")) %>%
    `[`(,c(ncol(.):1)) %>%
    as_tibble(.name_repair="unique") %>%
    transmute(method = ...3,
              expr = ...2,
              cohort = ...1,
              data = io$load(args$net)) %>%
    mutate(hits = purrr::map(data, check_hits, real=real)) %>%
    select(-data) %>%
    tidyr::unnest()

save(hits, file=args$outfile)

p = ggplot(hits, aes(x=FP, y=TP, color=method)) +
    geom_abline(data=ng, aes(slope=slope, intercept=0), color="grey", linetype="dashed") +
    geom_line() +
    facet_grid(expr ~ cohort) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
pdf(args$plotfile, 12, 8)
print(p)
dev.off()
