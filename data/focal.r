sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', 'TableS2D.xlsx'),
    opt('o', 'outfile', 'file to save to', 'focal.RData'))
