all:
	make tex && make bib && make tex
tex:
	pdflatex main
bib:
	biber main
wc:
	python3 ../wc/wc.py abstract s.intro s.prelim s.math s.appl s.disc
cls:
	cp ~/dev/tex/pkg/ms/ms/ms.cls .
clean:
	texclean .
