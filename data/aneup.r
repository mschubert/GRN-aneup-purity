sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', ''),
    opt('o', 'outfile', 'file to save to', 'aneup.RData'))
