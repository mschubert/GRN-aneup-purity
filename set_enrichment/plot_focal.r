library(dplyr)
library(cowplot)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'RData', 'merge.RData'),
    opt('r', 'region', 'character string', 'focal'),
    opt('p', 'plotfile', 'pdf', 'merge.pdf'))

allseg = filter(res, seg_id == "all", regions == args$region) %>% select(-seg_id)

pdf(args$plotfile, 12, 12)
ggplot(allseg, aes(x=edges, y=estimate, color=method)) +
    geom_hline(yintercept=1, color="black", linetype="dashed") +
    geom_line() +
    geom_point() +
    facet_grid(expr ~ cohort) +
    scale_x_log10() +
    scale_y_log10(breaks=c(0.1,1,3,10,30,100)) +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("Odds ratio for links in the same CNA")

ggplot(allseg, aes(x=edges, y=pmax(0,exp_fp), color=method)) +
    geom_line() +
    geom_point() +
    facet_grid(expr ~ cohort) +
    scale_x_log10() +
    scale_y_log10() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("Expected number of false positive links")

ggplot(allseg, aes(x=edges, y=fpr_total*100, color=method)) +
    geom_line() +
    geom_point() +
    facet_grid(expr ~ cohort, scales="free_y") +
    scale_x_log10() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("FPR genome-wide")

ggplot(allseg, aes(x=edges, y=pmax(0,fpr_seg*100), color=method)) +
    geom_line() +
    geom_point() +
    facet_grid(expr ~ cohort) +
    scale_x_log10() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("FPR within CNA")
dev.off()
