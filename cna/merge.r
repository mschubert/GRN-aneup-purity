library(dplyr)
library(cowplot)
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
    opt('p', 'plotfile', 'pdf', 'merge.pdf'),
    arg('infiles', 'FET files to load', arity='*',
        list.files(".", "\\.RData", recursive=TRUE, full.names=TRUE)))

res = do.call(rbind, strsplit(args$infiles, "[./_]")) %>%
    as.data.frame() %>%
    transmute(regions = V3,
              method = V4,
              cohort = V5,
              data = io$load(args$infiles)) %>%
    tidyr::unnest() %>%
    mutate(affected = purrr::map2(estimate, links, affected)) %>%
    select(-links) %>%
    tidyr::unnest()

pdf(args$plotfile)
ggplot(res) +
    geom_bar(aes(x=cohort, y=estimate), stat="identity") +
    facet_grid(method ~ regions, scales="free_y") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("Odds ratio for links in the same CNA")

ggplot(res) +
    geom_bar(aes(x=cohort, y=exp_fp), stat="identity") +
    facet_grid(method ~ regions, scales="free_y") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("Expected number of false positive links")

ggplot(res) +
    geom_bar(aes(x=cohort, y=fpr_total*100), stat="identity") +
    facet_grid(method ~ regions, scales="free") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("FPR genome-wide")

ggplot(res) +
    geom_bar(aes(x=cohort, y=fpr_seg*100), stat="identity") +
    facet_grid(method ~ regions, scales="free") +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("FPR within CNA")
dev.off()

save(res, file=args$outfile)