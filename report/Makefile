SRC = $(wildcard src/*.Rmd)

Supplements.pdf: Supplements.Rmd $(SRC)
	R --no-save --no-restore -e "rmarkdown::render('$<', 'pdf_document')"
