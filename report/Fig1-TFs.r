library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
io = import('io')

###
### Overview of genes, TFs, samples, etc. of DREAM4 vs. DREAM5 vs. TCGA
###
config = io$read_yaml("../config.yaml")
ng = io$load("../data/ng.RData")
tfa = io$load("../data/tf_annot.RData")
tfb = io$load("../data/tf_binding.RData")

dream = data.frame(
    collec = c(rep("DREAM 5", 3), rep("DREAM 4", 5)),
    net = c("In-silico (1)", "E. coli (2)", "S. cerevisiae (3)", paste("Network", 1:5)),
    ng = c(1643, 4511, 5950, rep(100, 5)),
    ntf = c(195, 334, 333, rep(100, 5)),
    nsmp = c(805, 805, 536, rep(100, 5)),
    nint = c(4012, 2066, 3940, 176, 249, 195, 211, 193)
)

lookup = c(
    ng = "# Genes",
    nsmp = "# Samples",
    ntf = "# TFs",
    nint = "# Interactions"
)

both = ng %>%
    transmute(collec = "TCGA",
              net = cohort,
              ng = ng, ntf=ntf, nsmp=nsmp, nint=nint) %>%
    bind_rows(dream) %>%
    tidyr::gather("field", "value", -net, -collec) %>%
    mutate(net = factor(net, levels=rev(c(config$cohorts, dream$net))),
           field = factor(lookup[field], levels=lookup))

p11 = ggplot(both, aes(y=net, x=value)) +
    ggstance::geom_colh() +
    geom_text(aes(label=value, x=0), hjust=-0.5) +
    facet_grid(collec~field, scales="free", space="free_y") +
    labs(tag = "b") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle=60, hjust=1))

###
### TF annot. vs. binding vs. all genes
###
ng = 22000
tfa = length(tfa)
tgb = length(unique(unlist(tfb)))
tfb = length(tfb)
rects = matrix(nrow=4, dimnames=list(c("xmin", "ymin", "xmax", "ymax"), NULL), c(
    0,0, ng,ng, # all gene conntections
    0,0+6000, ng,tfa+6000, # annotated TFs x genes
#    0,0, tgb,tfb, # binding TFs x bound genes
    0,0, 100, 100, # dream4
    0,0, dream[[1,"ng"]], dream[[1,"ntf"]], # dream5 in silico
    0,0, dream[[2,"ng"]], dream[[2,"ntf"]], # dream5 e. coli
    0,0, dream[[3,"ng"]], dream[[3,"ntf"]], # dream5 s. cerevisiae
NULL)) %>% t() %>% as.data.frame()
rects[3:6,c('xmin','xmax')] = rects[3:6,c('xmin','xmax')] + 1000
rects[3:6,c('ymin','ymax')] = rects[3:6,c('ymin','ymax')] + rep(1:4, length.out=8) * 1500 + 10000
rects = rects %>%
    mutate(color = letters[1:nrow(rects)],
           txt = c("All human genes", "Annotated human TFs x genes",
                   "DREAM 4", paste(dream$collec, dream$net)[1:3]),
           txtx = xmax + 1000,
           txty = mapply(function(x,y) mean(c(x,y)), ymin, ymax))
rects[1:2,'txtx'] = 1000
rects[1,'txty'] = 3000

p12 = ggplot(rects) +
    geom_rect(aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax, fill=color), color="black", alpha=0.8) +
    geom_text(aes(label=txt, x=txtx, y=txty), vjust=0.5, hjust=0) +
    coord_fixed() +
    guides(fill=FALSE) +
    labs(x = "Number of genes",
         y = "Number of regulators",
         tag = "a")

###
### Comparison of inferred networks with vs. without TFs
###
perf = io$load("../set_enrichment/TFbinding_enrichment.RData")
hits = filter(perf$hits, expr == "naive", size<=1e6) %>%
    mutate(mstr = factor(unlist(config$method_str[method]), levels=names(config$mcol)),
           has_tf = method %in% config$has_tf)
hl = data.frame(xmin=0, ymin=0, xmax=2e5, ymax=1200)
#hitsTF = filter(hits, method %in% config$has_tf, size<=hl$xmax)
hitsTF = filter(hits, size<=hl$xmax)
slope = filter(perf$ng, expr == "naive")

p2 = ggplot(hits, aes(x=size, y=TP, color=mstr)) +
    geom_rect(data=hl, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              inherit.aes=FALSE, fill="#F0E8D8", color="black", linetype="dotted") +
    geom_abline(data=slope, aes(slope=slope, intercept=0), color="grey", linetype="dashed") +
    geom_line(aes(size=has_tf, group=method), alpha=0.7) +
    facet_wrap(~cohort, nrow=1) +
    labs(x = "",
         y = "TF:TG recovered",
         tag = "c") +
    scale_size_manual(values=c(0.6,1)) +
    scale_color_manual(values=unlist(config$mcol), na.translate=FALSE) +
    guides(color = guide_legend(title="Method"),
           size = guide_legend(title="TF annotations")) +
    theme(axis.text.x = element_text(angle=45, hjust=1))

p3 = ggplot(hitsTF, aes(x=size, y=TP, color=mstr)) +
    geom_abline(data=slope, aes(slope=slope, intercept=0), color="grey", linetype="dashed") +
    geom_line(aes(size=has_tf, group=method), alpha=0.7) +
    ylim(c(0, hl$ymax)) +
    facet_wrap(~cohort, nrow=1) +
    labs(x = "# of edges considered",
         y = "TF:TG recovered",
         tag = "d") +
    scale_size_manual(values=c(0.6,1)) +
    scale_color_manual(values=unlist(config$mcol)) +
    guides(color=FALSE, size=FALSE) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
#          panel.background = element_rect(fill="#F0E8D8"))

pdf("Fig1-TFs.pdf", 14, 12)
({ (p12 | p11) + plot_layout(widths=c(1,2)) } / p2 / p3 ) + plot_layout(heights=c(3,2,2))
dev.off()
