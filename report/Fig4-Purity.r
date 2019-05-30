library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')

config = io$read_yaml("../config.yaml")

lookup = c(
    fpr_seg = "within purity genes",
    fpr_total = "fraction of total"
)

dset = io$load("../set_enrichment/merge.RData") %>%
    filter(expr == "naive", seg_id %in% c("purity1000", "purity5000")) %>%
    mutate(fpr_seg = pmax(0, fpr_seg)) %>%
    select(seg_id, method, cohort, edges, fpr_seg, fpr_total) %>%
    tidyr::gather("kind", "value", fpr_seg, fpr_total) %>%
    filter(! (kind == "fpr_seg" & method == "TFbinding"),
           method != "zscore.wrap") %>%
    mutate(kind = factor(lookup[kind], levels=lookup),
           mstr = unlist(config$method_str[method]),
           has_tf = method %in% config$has_tf)
p1000 = filter(dset, seg_id == "purity1000")
p5000 = filter(dset, seg_id == "purity5000")

p1 = ggplot(p1000, aes(x=edges, y=value*100, color=mstr)) +
    geom_line(aes(size=has_tf, group=method)) +
    geom_point(size=0.5) +
    facet_grid(kind ~ cohort, scales="free_y") +
    scale_x_log10() +
    scale_size_manual(values=c(0.6, 1)) +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    guides(color = guide_legend(title="Method"),
           size = guide_legend(title="TF annotations")) +
    labs(x = "# of edges",
         y = "% of edges expected FPs",
         tag = "a",
         title = "1000 genes most correlated with purity")

p2 = ggplot(p5000, aes(x=edges, y=value*100, color=mstr)) +
    geom_line(aes(size=has_tf, group=method)) +
    geom_point(size=0.5) +
    facet_grid(kind ~ cohort, scales="free_y") +
    scale_x_log10() +
    scale_size_manual(values=c(0.6,1)) +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    guides(color=FALSE, size=FALSE) +
    labs(x = "# of edges",
         y = "% of edges expected FPs",
         tag = "b",
         title = "5000 genes most correlated with purity")

pdf("Fig4-Purity.pdf", 14, 12)
p1 / p2
dev.off()
