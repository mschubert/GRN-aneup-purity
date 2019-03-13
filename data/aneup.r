io = import('io')
sys = import('sys')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', ''),
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('o', 'outfile', 'file to save to', 'aneup.RData'))

cohorts = io$read_yaml(args$config)$cohorts

chrs = lapply(cohorts, tcga$cna_chrs) %>%
    lapply(reshape2::melt, value.name="copies") %>%
    dplyr::bind_rows()

save(chrs, file=args$outfile)
