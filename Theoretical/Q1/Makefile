
OBJ = Q1.tex

doc: ${OBJ}
	@pdflatex $^
	@bibtex Q1.aux
	@pdflatex $^
	@pdflatex $^

.PHONY:clean
clean:
	rm -f *.aux *.log *.out *.synctex.gz  *.toc *.bbl *.run.xml *.blg 