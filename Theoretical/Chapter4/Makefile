OBJ = main.tex

all: main.pdf

main.pdf: ${OBJ}
	@xelatex $^
	@xelatex $^

doc: ${OBJ}
	@xelatex $^
	@xelatex $^

.PHONY:clean
clean:
	rm -f *.aux *.log *.out *.synctex.gz  *.toc *.bbl *.run.xml *.blg 