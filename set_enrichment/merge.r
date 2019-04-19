library(dplyr)
io = import('io')
sys = import('sys')

affected = function(estimate, lnk) { # also count genes affected?
    fpr_seg = 1 - 1/estimate
    exp_fp = fpr_seg * lnk["real", "seg"]
    fpr_total = exp_fp / (lnk["real", "all"] + lnk["real", "seg"]) # 'all' has FET margins
    data.frame(exp_fp=exp_fp, fpr_seg=fpr_seg, fpr_total=fpr_total)
}

args = sys$cmd$parse(
    opt('o', 'outfile', 'RData', 'merge.RData'),
    arg('infiles', 'FET files to load', arity='*',
        list.files(".", "\\.RData", recursive=TRUE, full.names=TRUE)))

res = do.call(rbind, strsplit(tools::file_path_sans_ext(args$infiles), "[/_]")) %>%
    as.data.frame() %>%
    transmute(regions = V2,
              method = V3,
              expr = V4,
              cohort = V5,
              data = io$load(args$infiles)) %>%
    tidyr::unnest() %>%
    mutate(affected = purrr::map2(estimate, links, affected)) %>%
    select(-links) %>%
    tidyr::unnest()

save(res, file=args$outfile)
