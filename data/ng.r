library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('a', 'annot', 'TF annotation RData', 'tf_annot.RData'),
    opt('o', 'outfile', 'file to save df to', 'ng.RData'),
    arg('expr', 'expr RData', arity='*',
        list.files("expr_naive", recursive=TRUE, full.names=TRUE)))

#TODO: filter psbl_tf by actually expressed?
ntf = length(io$load(args$annot))

ng = do.call(rbind, strsplit(tools::file_path_sans_ext(args$expr), "[_/]")) %>%
    `[`(,c(ncol(.):1)) %>%
    as_tibble(.name_repair="unique") %>%
    transmute(cohort = ...1,
              expr = ...2,
              ng = sapply(args$expr, function(e) nrow(io$load(e))),
              psbl = 0.5 * (ng-1) * ng,
              psbl_tf = ntf * ng - 0.5 * (ntf-1) * ntf)

save(ng, file=args$outfile)
