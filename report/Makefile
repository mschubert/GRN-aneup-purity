R = R --no-save --no-restore -e
SRC = $(wildcard *.r)
PDF = $(SRC:%.r=%.pdf)

all: $(PDF) $(PDF:%.pdf=%.png) Figures.pdf

Figures.pdf:
	pdftk $(sort $(PDF)) cat output $@

%.pdf: %.r
	Rscript $<

%.pdf: %.Rmd
	$(R) "rmarkdown::render('$<', 'pdf_document')"

%.png: %.pdf
	convert -density 400 $< -resize 25% -flatten $@
