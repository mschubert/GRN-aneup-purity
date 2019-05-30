library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')

config = io$read_yaml("../config.yaml")

lookup = c(
    fpr_seg = "within CNA",
    fpr_total = "fraction of total"
)

dset = io$load("../set_enrichment/merge.RData") %>%
    filter(expr == "naive", seg_id == "all") %>%
    mutate(fpr_seg = pmax(0, fpr_seg)) %>%
    select(regions, method, cohort, edges, fpr_seg, fpr_total) %>%
    tidyr::gather("kind", "value", fpr_seg, fpr_total) %>%
    filter(! (kind == "fpr_seg" & method == "TFbinding"),
           method != "zscore.wrap") %>%
    mutate(kind = factor(lookup[kind], levels=lookup),
           mstr = unlist(config$method_str[method]),
           has_tf = method %in% config$has_tf)
focal = filter(dset, regions == "focal")
aneup = filter(dset, regions == "aneup")

levels(aneup$kind)[1] = "within Chromosome"

p1 = ggplot(focal, aes(x=edges, y=value*100, color=mstr)) +
    geom_line(aes(size=has_tf, group=method), alpha=0.7) +
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
         title = "Focal CNAs (as defined by RACS in GDSC)")

p2 = ggplot(aneup, aes(x=edges, y=value*100, color=mstr)) +
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
         title = "Aneuploidies (chromosome copy number changes)")

pdf("Fig3-CNAs.pdf", 14, 12)
p1 / p2
dev.off()
