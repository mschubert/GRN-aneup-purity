library(dplyr)
library(magrittr)
io = import('io')
sys = import('sys')
seq = import('seq')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', ''),
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('o', 'outfile', 'file to save to', 'aneup.RData'))

cohorts = io$read_yaml(args$config)$cohorts

copies = lapply(cohorts, tcga$cna_chrs) %>%
    setNames(cohorts)

genes = seq$coords$gene("hgnc_symbol", chromosomes=c(1:22,'X')) %>%
    select(chromosome_name, gene=external_gene_name) %>%
    group_by(chromosome_name) %>%
    tidyr::nest() %>%
    arrange(chromosome_name) %$%
    setNames(data, chromosome_name)

save(copies, genes, file=args$outfile)
