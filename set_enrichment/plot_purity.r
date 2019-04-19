library(dplyr)
io = import('io')
sys = import('sys')
plt = import('./plot')
eplot = import('plot')

args = sys$cmd$parse(
    opt('i', 'infile', 'RData', 'merge.RData'),
    opt('r', 'region', 'character string', 'purity'),
    opt('p', 'plotfile', 'pdf', 'plot_purity.pdf'))

allseg = io$load(args$infile) %>%
    filter(regions == args$region) %>%
    split(.$seg_id)

pdf(args$plotfile, 12, 8)
for (segs in allseg) {
    eplot$error(segs$seg_id[1])
    print(plt$odds(segs))
    print(plt$fp_total(segs))
    print(plt$fpr_genome(segs))
    print(plt$fpr_segment(segs))
}
dev.off()
