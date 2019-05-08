library(dplyr)
io = import('io')
sys = import('sys')

count_cohort = function(fname) {
    expr = io$load(fname)
    data.frame(ng = nrow(expr),
         ntf = length(intersect(rownames(expr), tfa)),
         nsmp = ncol(expr),
         nint = nrow(filter(tftg, ind %in% rownames(expr) & values %in% rownames(expr)))
    )
}

args = sys$cmd$parse(
    opt('a', 'annot', 'TF annotation RData', 'tf_annot.RData'),
    opt('b', 'binding', 'TF binding RData', 'tf_binding.RData'),
    opt('o', 'outfile', 'file to save df to', 'ng.RData'),
    arg('expr', 'expr RData', arity='*',
        list.files("expr_naive", recursive=TRUE, full.names=TRUE)))

tfa = io$load(args$annot)
tfb = io$load(args$binding)
tftg = rev(stack(tfb))

ng = do.call(rbind, strsplit(tools::file_path_sans_ext(args$expr), "[_/]")) %>%
    `[`(,c(ncol(.):1)) %>%
    as_tibble(.name_repair="unique") %>%
    transmute(cohort = ...1,
              expr = ...2,
              values = purrr::map(args$expr, count_cohort)) %>%
    tidyr::unnest() %>%
    mutate(psbl = 0.5 * (ng-1) * ng,
           psbl_tf = ntf * ng - 0.5 * (ntf-1) * ntf)

save(ng, file=args$outfile)
