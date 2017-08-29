syn:	syn.tex
	pdflatex syn
	bibtex syn
	pdflatex syn
	pdflatex syn
	open syn.pdf
figs:
	rm -f figs.eps
	latex figs
	dvips -o figs.ps figs.dvi
	ps2eps figs.ps
	ps2pdf figs.ps
	open figs.pdf
clean:
	rm -f *.bbl *.blg syn.pdf *.dvi *.log *.aux *Notes.bib *.ps *.nav *.snm *.toc *.out
