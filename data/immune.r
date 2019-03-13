sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', '1-s2.0-S1074761318301213-mmc2.xlsx'),
    opt('o', 'outfile', 'file to save to', 'immune.RData'))
