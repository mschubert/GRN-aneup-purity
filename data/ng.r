library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'outfile', 'file to save df to', 'ng.RData'),
    arg('expr', 'expr RData', arity='*',
        list.files("expr_naive", recursive=TRUE, full.names=TRUE)))

ng = do.call(rbind, strsplit(tools::file_path_sans_ext(args$expr), "[_/]")) %>%
    `[`(,c(ncol(.):1)) %>%
    as_tibble(.name_repair="unique") %>%
    transmute(cohort = ...1,
              expr = ...2,
              ng = sapply(args$expr, function(e) nrow(io$load(e))),
              psbl = 0.5 * (ng-1) * ng)

save(ng, file=args$outfile)
