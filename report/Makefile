R = R --no-save --no-restore -e
SRC = $(wildcard *.r)

all: $(SRC:%.r=%.pdf)

%.pdf: %.r
	Rscript $<

%.pdf: %.Rmd
	$(R) "rmarkdown::render('$<', 'pdf_document')"
