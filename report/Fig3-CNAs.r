library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')

lookup = c(
    fpr_seg = "FPR within CNA",
    fpr_total = "FPR over all edges"
)

dset = io$load("../set_enrichment/merge.RData") %>%
    filter(expr == "naive", seg_id == "all") %>%
    mutate(fpr_seg = pmax(0, fpr_seg)) %>%
    select(regions, method, cohort, edges, fpr_seg, fpr_total) %>%
    tidyr::gather("kind", "value", fpr_seg, fpr_total) %>%
    filter(! (kind == "fpr_seg" & method == "TFbinding"),
           method != "zscore.wrap") %>%
    mutate(kind = factor(lookup[kind], levels=lookup))
focal = filter(dset, regions == "focal")
aneup = filter(dset, regions == "aneup")

levels(aneup$kind)[1] = "FPR within Chromosome"

p1 = ggplot(focal, aes(x=edges, y=value*100, color=method)) +
    geom_line() +
    geom_point(size=0.5) +
    facet_grid(kind ~ cohort, scales="free_y") +
    scale_x_log10() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x = "# of edges",
         y = "% of edges",
         tag = "a",
         title = "Focal CNAs (as defined by RACS in GDSC)")

p2 = ggplot(aneup, aes(x=edges, y=value*100, color=method)) +
    geom_line() +
    geom_point(size=0.5) +
    facet_grid(kind ~ cohort, scales="free_y") +
    scale_x_log10() +
    theme(panel.grid.major.y = element_line(color="grey", linetype="dashed"),
          axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x = "# of edges",
         y = "% of edges",
         tag = "b",
         title = "Aneuploidies (chromosome copy number changes)")

pdf("Fig3-CNAs.pdf", 14, 12)
p1 / p2
dev.off()
