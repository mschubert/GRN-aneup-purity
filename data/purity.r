sys = import('sys')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', ''),
    opt('o', 'outfile', 'file to save to', 'purity.RData'))

purity = tcga$purity()
save(purity, file=args$outfile)
