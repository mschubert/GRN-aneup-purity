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

allseg = filter(res, seg_id == "all") %>% select(-seg_id)

pdf(args$plotfile, 12, 12)
ggplot(allseg, aes(x=edges, y=estimate, color=method)) +
    geom_hline(yintercept=1, color="black", linetype="dashed") +
    geom_line() +
    geom_point() +
    facet_grid(regions+expr ~ cohort) +
    scale_x_log10() +
    scale_y_log10(breaks=c(0.1,1,3,10,30,100)) +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("Odds ratio for links in the same CNA")

ggplot(allseg, aes(x=edges, y=exp_fp, color=method)) +
    geom_line() +
    geom_point() +
    facet_grid(regions+expr ~ cohort) +
    scale_x_log10() +
    scale_y_log10() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("Expected number of false positive links")

ggplot(allseg, aes(x=edges, y=fpr_total*100, color=method)) +
    geom_line() +
    geom_point() +
    facet_grid(regions+expr ~ cohort, scales="free_y") +
    scale_x_log10() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("FPR genome-wide")

ggplot(allseg, aes(x=edges, y=pmax(0,fpr_seg*100), color=method)) +
    geom_line() +
    geom_point() +
    facet_grid(regions+expr ~ cohort) +
    scale_x_log10() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("FPR within CNA")
dev.off()

save(res, file=args$outfile)
