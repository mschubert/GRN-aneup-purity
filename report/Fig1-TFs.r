library(dplyr)
library(cowplot)
io = import('io')

ng = io$load("../data/ng.RData")
tfa = io$load("../data/tf_annot.RData")
tfb = io$load("../data/tf_binding.RData")

dream = data.frame(
    collec = c(rep("DREAM 5", 3), rep("DREAM 4", 5)),
    net = c("In-silico (1)", "E. coli (2)", "S. cervisiae (3)", paste("Network", 1:5)),
    ng = c(1643, 4511, 5950, rep(100, 5)),
    ntf = c(195, 334, 333, rep(100, 5)),
    nsmp = c(805, 805, 536, rep(100, 5)),
    nint = c(4012, 2066, 3940, 176, 249, 195, 211, 193)
)

both = ng %>%
    transmute(collec = "TCGA",
              net = cohort,
              ng = ng, ntf=ntf, nsmp=nsmp, nint=nint) %>%
    bind_rows(dream) %>%
    tidyr::gather("field", "value", -net, -collec)

p = ggplot(both, aes(y=net, x=value)) +
    ggstance::geom_colh() +
    geom_text(aes(label=value, x=0), hjust=-0.5) +
    facet_grid(collec~field, scales="free", space="free_y")

pdf("Fig1-TFs.pdf", 12, 6)
print(p)
dev.off()
