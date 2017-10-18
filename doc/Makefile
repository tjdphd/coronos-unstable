pdf:	syn.tex
	pdflatex syn
	bibtex syn
	pdflatex syn
	pdflatex syn
	open syn.pdf
t4:	syn.tex
	mk4ht htlatex syn
	bibtex syn
	mk4ht htlatex syn
figs:
	rm -f figs.eps
	latex figs
	dvips -o figs.ps figs.dvi
	ps2eps figs.ps
	ps2pdf figs.ps
	open figs.pdf
clean:
	rm -f *.bbl *.blg syn.pdf *.dvi *.log *.aux *Notes.bib *.ps *.nav *.snm *.toc *.out	\
	syn.aux syn.log syn.out syn.pdf								\
	syn.toc syn.bbl syn.blg syn.dvi								\
	syn.idv syn.lg syn.epub syn.html							\
	syn.4ct syn.4tc syn.css syn.tmp								\
	syn.xref syn*.png texput.log t4.html							\
	syn.glg syn.glo ebook.gls  syn.ist							\
	syn.tlg
