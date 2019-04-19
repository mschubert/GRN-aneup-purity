library(dplyr)
io = import('io')
sys = import('sys')
plt = import('./plot')

args = sys$cmd$parse(
    opt('i', 'infile', 'RData', 'merge.RData'),
    opt('r', 'region', 'character string', 'focal'),
    opt('p', 'plotfile', 'pdf', 'merge.pdf'))

allseg = io$load(args$infile) %>%
    filter(seg_id == "all", regions == args$region) %>% select(-seg_id)

pdf(args$plotfile, 12, 8)
plt$odds(allseg)
plt$fp_total(allseg)
plt$fpr_genome(allseg)
plt$fpr_segment(allseg)
dev.off()
