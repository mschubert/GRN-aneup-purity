library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')

lookup = c(
    fpr_seg = "FPR within purity genes",
    fpr_total = "FPR over all edges"
)

dset = io$load("../set_enrichment/merge.RData") %>%
    filter(expr == "naive", seg_id %in% c("purity500", "purity3000")) %>%
    mutate(fpr_seg = pmax(0, fpr_seg)) %>%
    select(seg_id, method, cohort, edges, fpr_seg, fpr_total) %>%
    tidyr::gather("kind", "value", fpr_seg, fpr_total) %>%
    filter(! (kind == "fpr_seg" & method == "TFbinding"),
           method != "zscore.wrap") %>%
    mutate(kind = factor(lookup[kind], levels=lookup))
p500 = filter(dset, seg_id == "purity500")
p3000 = filter(dset, seg_id == "purity3000")

p1 = ggplot(p500, aes(x=edges, y=value*100, color=method)) +
    geom_line() +
    geom_point(size=0.5) +
    facet_grid(kind ~ cohort, scales="free_y") +
    scale_x_log10() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x = "# of edges",
         y = "% of edges",
         tag = "a",
         title = "500 genes most correlated with purity")

p2 = ggplot(p3000, aes(x=edges, y=value*100, color=method)) +
    geom_line() +
    geom_point(size=0.5) +
    facet_grid(kind ~ cohort, scales="free_y") +
    scale_x_log10() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x = "# of edges",
         y = "% of edges",
         tag = "b",
         title = "3000 genes most correlated with purity")

pdf("Fig4-Purity.pdf", 14, 12)
p1 / p2
dev.off()
