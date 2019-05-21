library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')
tcga = import('data/tcga')

config = io$read_yaml("../config.yaml")
focal = io$load("../data/focal.RData")
aneup = io$load("../data/aneup.RData")
purity = io$load('../data/purity.RData')
hl = config$cohorts

###
### Top row overview across all cohorts
###
pur = lapply(purity$estimate, reshape2::melt) %>%
    bind_rows(.id = "cohort") %>%
    transmute(cohort = cohort,
              Sample = Var1,
              purity = value,
              hl = factor(cohort, levels=hl))
cohorts = unique(pur$cohort)
anp = aneup$samples %>%
    filter(cohort %in% cohorts,
           substr(Sample, 14, 15) %in% c("01", "02", "03")) %>% # primary, recurrent, metastasis
    mutate(hl = factor(cohort, levels=hl))
#    tcga$filter(along="Sample", cancer=TRUE)
focs = focal$sums %>%
    filter(cohort %in% cohorts | cohort == "COAD/READ") %>%
    mutate(pct = 100 * bases / 3.234e9,
           hl = factor(cohort, levels=hl))

p11 = ggplot(focs, aes(x=n_regions, y=pct)) +
    geom_point(aes(color=hl), size=2) +
    ggrepel::geom_text_repel(aes(label=cohort)) +
    guides(color=FALSE) +
    labs(x = "Number of regions",
         y = "% genome covered",
         tag = "a")
p12 = ggplot(anp, aes(x=forcats::fct_reorder(cohort, aneuploidy), y=aneuploidy)) +
    geom_boxplot(aes(fill=hl), outlier.shape=NA) +
    guides(fill=FALSE) +
    labs(y = "Aneuploidy",
         tag = "b") +
    theme(axis.text.x = element_text(angle=90, hjust=1),
                 axis.title.x = element_blank())
p13 = ggplot(pur, aes(x=forcats::fct_reorder(cohort, purity), y=purity)) +
    geom_bar(aes(fill=hl), stat="summary", fun.y="mean") +
    stat_summary(fun.data=function(x) mean_sdl(x, mult=1), geom="errorbar", width=0.5) +
    ylim(0, 1) +
    labs(y = "Purity",
         tag = "c") +
    guides(fill=FALSE) +
    theme(axis.text.x = element_text(angle=90, hjust=1),
                 axis.title.x = element_blank())

###
### Focal sample distribution
###
df = lapply(focal$estimate[config$cohorts], reshape2::melt) %>%
    bind_rows(.id="cohort") %>%
    mutate(name = sub("[A-Z]+", "", name),
           name = factor(name, levels=gtools::mixedsort(unique(name))))

ng = lapply(focal$sets[config$cohorts], function(x) sapply(x, length)) %>%
    lapply(stack) %>%
    bind_rows(.id = "cohort") %>%
    transmute(cohort=cohort, ng=values, name=sub("[A-Z]+", "", ind)) %>%
    filter(paste(cohort, name) %in% paste(df$cohort, df$name))

y_labels = seq_along(levels(df$name))
y_labels[(y_labels-1) %% 5 != 0] = ""

p21 = ggplot(df) +
    geom_vline(xintercept=1:9, linetype="dotted", color="grey") +
    ggstance::geom_boxploth(aes(x=value, y=name, fill=cohort),
                            outlier.shape=1, outlier.size=0.2, outlier.alpha=0.5) +
    facet_grid(cohort ~ ., scales="free", space="free_y") +
    guides(fill=FALSE) +
    scale_x_continuous(trans="log2", breaks=1:9) +
    scale_y_discrete(labels=y_labels) +
    labs(tag = "d",
         x = "DNA copies of segment",
         y = "Regions of focal amplification/deletion (RACS)")

p22 = ggplot(ng) +
    ggstance::geom_barh(aes(x=ng, y=name, fill=cohort), color="darkgrey", stat="identity") +
    facet_grid(cohort ~ ., scales="free", space="free_y") +
    guides(fill=FALSE) +
    scale_x_continuous(trans="log2") +
    scale_y_discrete(labels=y_labels) +
    labs(x = "# genes") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())

p2 = p21 + p22 + plot_layout(nrow=1, widths=c(5, 1)) & theme(plot.margin = unit(c(0,0,0,0), "cm"))

###
### Aneuploidy sample distribution
###
df = lapply(aneup$estimate[config$cohorts], reshape2::melt) %>%
    bind_rows(.id="cohort") %>%
    mutate(chrs = factor(chrs, levels=c(1:22,'X')))

ng = sapply(aneup$sets, length) %>%
    stack() %>%
    transmute(ng=values, chrs=ind) %>%
    inner_join(distinct(df[c("cohort", "chrs")]))

y_labels = seq_along(levels(df$chrs))
y_labels[(y_labels-1) %% 5 != 0 | y_labels > 20] = ""
y_labels[length(y_labels)] = "X"

p31 = ggplot(df) +
    geom_vline(xintercept=2:3, linetype="dotted", color="grey") +
    ggstance::geom_boxploth(aes(x=value, y=chrs, fill=cohort),
                            outlier.shape=1, outlier.size=0.2, outlier.alpha=0.5) +
    facet_grid(cohort ~ ., scales="free", space="free_y") +
    guides(fill=FALSE) +
    scale_x_continuous(trans="log2", breaks=2:3) +
    scale_y_discrete(labels=y_labels) +
    labs(tag = "e",
         x = "Ploidy (number of chromosomes)",
         y = "Chromosome")

p32 = ggplot(ng) +
    ggstance::geom_barh(aes(x=ng, y=chrs, fill=cohort), color="darkgrey", stat="identity") +
    facet_grid(cohort ~ ., scales="free", space="free_y") +
    guides(fill=FALSE) +
    scale_y_discrete(labels=y_labels) +
    labs(x = "# genes") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())

p3 = p31 + p32 + plot_layout(nrow=1, widths=c(5, 1)) & theme(plot.margin = unit(c(0,0,0,0), "cm"))

###
### Assemble
###
pdf("Fig2-Distributions.pdf", 12, 16)
(p11 | p12 | p13) /
(p2 | p3) + plot_layout(heights=c(1,3))
dev.off()
