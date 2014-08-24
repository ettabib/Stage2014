.PHONY: all, update

CC = 
CFLAGS =  
LDFLAGS = 
EXEC =  
SRC=$(wildcard *.tex)
OBJ = $(SRC:.tex=.pdf) 

all: update

debug:
	@echo "OBJ=$(OBJ)"

update:rapport.tex
	@-pdflatex -interaction=nonstopmode rapport.tex > /dev/null 2> /dev/null
	@-cecho -yellow "-------< First compilation done >-------"
	@echo " "
	@-pdflatex -interaction=nonstopmode rapport.tex > /dev/null 2> /dev/null
	@-cecho -yellow "-------< Second compilation done >-------"
	@-echo " "
	@-makeindex rapport.idx
	@cecho -yellow "-------< Index compilation done >-------"
	@-pdflatex -interaction=nonstopmode rapport.tex > /dev/null 2> /dev/null
	@-cecho -yellow "-------< Third compilation done >-------"
	@-echo " "
	@-bibtex biblio.bib > /dev/null 2> /dev/null
	@-cecho -green "-------< Compiling bibtex done >-------"	
	@-echo " "
	@date | awk '{print $5}'


clean:
	@cecho  "***** Latex Clean *****"
	@-rm -f rapport.pdf > /dev/null 2> /dev/null
	@-rm -rf *.bbl *.idx *.blg *.toc *.log *.rel *.fmt *.prv *.aux > /dev/null 2> /dev/null

view:rapport.pdf
	@open rapport.pdf

mrproper: clean
