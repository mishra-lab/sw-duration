all:
	make tex && make bib && make tex
tex:
	pdflatex main
bib:
	biber main
wc:
	python3 ../wc/wc.py abstract intro pre math app discuss
cls:
	cp ~/dev/tex/pkg/ms/ms/ms.cls .
clean:
	texclean .
