io = import('io')
sys = import('sys')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', ''),
    opt('o', 'outfile', 'file to save to', 'purity.RData'))

purity = tcga$purity() %>%
    filter(!is.na(estimate))

estimate = matrix(purity$estimate, ncol=1, dimnames=list(purity$Sample, "purity")) %>%
    narray::split(along=1, subsets=purity$cohort)

save(purity, file=args$outfile)
