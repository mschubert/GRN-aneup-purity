library(dplyr)
library(cowplot)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'outfile', 'RData', 'merge.RData'),
    opt('p', 'plotfile', 'pdf', 'merge.pdf'),
    arg('infiles', 'FET files to load', arity='*',
        list.files(".", "\\.RData", recursive=TRUE, full.names=TRUE)))

files = do.call(rbind, strsplit(args$infiles, "[./_]")) %>%
    as.data.frame() %>%
    transmute(regions = V3,
              method = V4,
              cohort = V5,
              data = io$load(args$infiles)) %>%
    tidyr::unnest()

pdf(args$plotfile)
ggplot(files) +
    geom_bar(aes(x=cohort, y=estimate), stat="identity") +
    facet_grid(method ~ regions, scales="free_y") +
    scale_y_log10() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    ggtitle("Odds ratio for links in the same CNA")
dev.off()

save(files, file=args$outfile)
