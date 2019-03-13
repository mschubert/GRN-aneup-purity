sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'outfile', 'file to save to', 'meta.RData'),
    opt('c', 'config', 'yaml', '../config.yaml'),
    arg('infiles', 'files to merge', arity='*',
        c("focal.RData", "aneup.RData", "purity.RData", "immune.RData")))
