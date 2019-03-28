library(dplyr)
io = import('io')
sys = import('sys')
tcga = import('data/tcga')

cohort_cnas = function(cohort, segs) {
    if (cohort %in% c("COAD", "READ"))
        cohort2 = "COAD/READ"
    else
        cohort2 = cohort

    message(cohort, " (segments ", cohort2, ")")

    cseg = segs %>%
        filter(Cancer.Type == cohort2) %>%
        select(chr, start, stop, name=Identifier) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)

    tcga$cna_custom(cohort, "name", cseg)
}

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', 'TableS2D.xlsx'),
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('o', 'outfile', 'file to save to', 'focal.RData'))

cohorts = io$read_yaml(args$config)$cohorts

segs = readxl::read_xlsx(args$infile, skip=19)
colnames(segs) = make.names(colnames(segs))

copies = lapply(cohorts, cohort_cnas, segs=segs) %>%
    setNames(cohorts)
genes = segs$Contained.genes %>%
    strsplit(",", fixed=TRUE) %>%
    setNames(segs$Identifier)

save(copies, genes, file=args$outfile)
