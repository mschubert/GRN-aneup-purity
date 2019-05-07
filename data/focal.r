library(dplyr)
io = import('io')
sys = import('sys')
tcga = import('data/tcga')

cohort_cnas = function(cohort, segs) {
    message(cohort)
    cseg = GenomicRanges::makeGRangesFromDataFrame(segs, keep.extra.columns=TRUE)
    tcga$cna_custom(cohort, "name", cseg)
}

cohort_genes = function(segs) {
    sets = strsplit(segs$Contained.genes, ",", fixed=TRUE) %>%
        setNames(segs$name)
}

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', 'TableS2D.xlsx'),
    opt('o', 'outfile', 'file to save to', 'focal.RData'))

segs = readxl::read_xlsx(args$infile, skip=19) %>%
    setNames(make.names(colnames(.))) %>%
    select(name=Identifier, cohort=Cancer.Type, chr, start, stop, Contained.genes)

sums = segs %>%
    group_by(cohort) %>%
    summarize(n_regions = dplyr::n(),
              bases = sum(abs(start - stop)))

segs = split(segs, segs$cohort)
segs$COAD = segs$READ = segs$`COAD/READ`
segs$`COAD/READ` = NULL
segs$PANCAN = NULL # can not get CNVs using TCGA API here

estimate = mapply(cohort_cnas, cohort=names(segs), segs=segs)
sets = lapply(segs, cohort_genes)

save(sums, estimate, sets, file=args$outfile)
