
OBJ = $(wildcard *.tex)

doc: ${OBJ}
	@pdflatex $^
	@bibtex main.aux
	@pdflatex $^
	@pdflatex $^

.PHONY:clean
clean:
	rm -f *.aux *.log *.out *.synctex.gz  *.toc *.bbl *.run.xml *.blg 