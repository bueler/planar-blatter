ALL: blatter.pdf

#### for PDF notes:

figures := profile.pdf bdryblatter.pdf conform.pdf sparsepattern.pdf q1mesh.pdf chapeau.png showtests.pdf
figures := $(addprefix figs/,$(figures))

blatter.pdf: blatter.tex ${figures}
	pdflatex blatter.tex
	bibtex blatter
	pdflatex blatter.tex
	pdflatex blatter.tex

#### utilities

.PHONY : clean

clean :
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *~

