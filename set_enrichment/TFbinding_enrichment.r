library(dplyr)
library(magrittr)
library(cowplot)
io = import('io')
sys = import('sys')
st = import('stats')

check_hits = function(net, real) {
    re = c(paste(net[[1]], net[[2]]), paste(net[[2]], net[[1]])) %in%
        paste(real[[1]], real[[2]]) + 0
    roc = data.frame(TP = cumsum(re),
                     FP = cumsum(1-re))
    roc2 = roc[re != 0 | c(re[-1], 1) != 0,]
}

args = sys$cmd$parse(
    opt('r', 'real', 'RData', '../networks/TFbinding/naive/ACC.RData'),
    opt('o', 'outfile', '.RData', 'TFbinding_background.RData'),
    opt('p', 'plotfile', 'pdf', 'TFbinding_background.pdf'),
    arg('net', 'networks', arity='*', '../networks/aracne/naive/ACC.RData'))

real = unique(io$load(args$real))

hits = do.call(rbind, strsplit(tools::file_path_sans_ext(args$net), "/")) %>%
    as_tibble(.name_repair="unique") %>%
    transmute(method = ...3,
              expr = ...4,
              cohort = ...5,
              data = io$load(args$net, drop=FALSE)) %>%
    mutate(hits = purrr::map(data, check_hits, real=real)) %>%
    select(-data) %>%
    tidyr::unnest()

save(hits, file=args$outfile)

pdf(args$plotfile, 12, 8)
ggplot(hits, aes(x=FP, y=TP, color=method)) +
    geom_line() +
    facet_grid(expr ~ cohort)
dev.off()
