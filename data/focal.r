library(dplyr)
io = import('io')
sys = import('sys')
tcga = import('data/tcga')

cohort_cnas = function(cohort, segs) {
    cseg = segs %>%
        filter(Cancer.Type == cohort) %>%
        select(chr, start, stop, name=Identifier, genes=Contained.genes) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)

    rng = tcga$cna_custom(cohort, "name", cseg, as_array=FALSE) %>%
        mutate(gene = strsplit(cseg$genes, ",", fixed=TRUE))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', 'TableS2D.xlsx'),
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('o', 'outfile', 'file to save to', 'focal.RData'))

cohorts = io$read_yaml(args$config)$cohorts

segs = readxl::read_xlsx(args$infile, skip=19)
colnames(segs) = make.names(colnames(segs))

copies = lapply(cohorts, cohort_cnas, segs=segs) %>%
    dplyr::bind_rows()
save(copies, file=args$outfile)
