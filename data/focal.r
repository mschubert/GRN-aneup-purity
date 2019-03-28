library(dplyr)
io = import('io')
sys = import('sys')
tcga = import('data/tcga')

cohort_cnas = function(segs) {
    cseg = GenomicRanges::makeGRangesFromDataFrame(segs, keep.extra.columns=TRUE)
    tcga$cna_custom(unique(segs$cohort), "name", cseg)
}

cohort_genes = function(segs) {
    sets = strsplit(segs$Contained.genes, ",", fixed=TRUE) %>%
        setNames(segs$name)
}

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', 'TableS2D.xlsx'),
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('o', 'outfile', 'file to save to', 'focal.RData'))

cohorts = io$read_yaml(args$config)$cohorts

segs = readxl::read_xlsx(args$infile, skip=19) %>%
    setNames(make.names(colnames(.))) %>%
    select(name=Identifier, cohort=Cancer.Type, chr, start, stop, Contained.genes) %>%
    split(.$cohort)
segs$COAD = segs$READ = segs$`COAD/READ`
segs$`COAD/READ` = NULL

estimate = lapply(segs, cohort_cnas)
sets = lapply(segs, cohort_genes)

save(estimate, sets, file=args$outfile)
