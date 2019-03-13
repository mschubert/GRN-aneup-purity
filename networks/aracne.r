sys = import('sys')
aracne = import('tools/aracne')

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('m', 'method', 'method identifier', 'aracne'),
    opt('o', 'outfile', '.RData to save to', 'aracne/ACC.RData'))

# infer nethods using different methods here
# save them as data.frames or graph objects
