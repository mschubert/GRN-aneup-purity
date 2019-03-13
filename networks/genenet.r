sys = import('sys')
gnet = import('tools/genenet')

args = sys$cmd$parse(
    opt('c', 'cohort', 'cohort identifier', 'ACC'),
    opt('m', 'method', 'method identifier', 'genenet'),
    opt('o', 'outfile', '.RData to save to', 'genenet/ACC.RData'))

# infer nethods using different methods here
# save them as data.frames or graph objects
