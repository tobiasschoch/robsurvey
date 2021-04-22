REM Pushd C:\My\useR\univariate
Rscript -e "Sweave('C:/My/code/robsurvey/inst/comparison/comparison.Rnw')"
REM Pushd C:\My\useR
pdflatex.exe .\comparison.tex


