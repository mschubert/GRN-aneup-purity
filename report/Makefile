R = R --no-save --no-restore -e
SRC = $(wildcard *.r)

all: $(SRC:%.r=%.pdf) $(SRC:%.r=%.png)

%.pdf: %.r
	Rscript $<

%.pdf: %.Rmd
	$(R) "rmarkdown::render('$<', 'pdf_document')"

%.png: %.pdf
	convert -density 400 $< -resize 25% -flatten $@
