cat $1.tex > which.aux
pdflatex main ;\
cp main.pdf ../../out/fig/tikz/$1.pdf
