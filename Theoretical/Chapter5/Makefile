OBJ = HW5_3220104819_WangHao.tex

all: HW5_3220104819_WangHao.pdf

HW5_3220104819_WangHao.pdf: ${OBJ}
	@xelatex $^
	@xelatex $^

doc: ${OBJ}
	@xelatex $^
	@xelatex $^

.PHONY:clean
clean:
	rm -f *.aux *.log *.out *.synctex.gz  *.toc *.bbl *.run.xml *.blg 