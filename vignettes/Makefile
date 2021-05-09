DOCUMENT=robsurvey
RSCRIPT=Rscript.exe

$(DOCUMENT).html:
	$(RSCRIPT) -e "library(rmarkdown); rmarkdown::render('$(DOCUMENT).Rmd', 'html_document')"

.PHONY: clean
clean:
	rm *.html
