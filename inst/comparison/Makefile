DOCUMENT=comparison
RSCRIPT=Rscript.exe
PDFLATEX=pdflatex.exe
PDFVIEWER=SumatraPDF.exe
R_LIB='C:/Rpackages/win-library/4.0'

all: sweave pdf1 pdf2

# pdflatex
pdf1: $(DOCUMENT).tex
	$(PDFLATEX) $(DOCUMENT).tex

pdf2: $(DOCUMENT).tex
	$(PDFLATEX) $(DOCUMENT).tex

# Sweave
sweave: $(DOCUMENT).Rnw
	$(RSCRIPT) -e ".libPaths($(R_LIB))" -e "Sweave('$(DOCUMENT).Rnw')"

# view pdf
view: $(DOCUMENT).pdf
	$(PDFVIEWER) $(DOCUMENT).pdf &

.PHONY:
final:
	rm *.aux *.out *.log *.tex

.PHONY:
clean:
	rm *.pdf *.aux *.out *.log *.tex
