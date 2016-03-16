all: paper

paper:
	pdflatex main
	bibtex main
	pdflatex main

clean:
	rm -f main.aux
	rm -f main.bbl
	rm -f main.blg
	rm -f main.log
	rm -f main.pdf
	rm -f main.out
	rm -f texput.log
