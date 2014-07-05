.PHONY: all


all:  update

debug:
	@echo "DIR_LATEX=$(DIR_LATEX)"

clean:
	@cecho --orange "***** Cleaning latex *****"
	@cd $(DIR_LATEX) && make clean

view:
	ln -s latex/rapport.pdf .
	@open rapport.pdf

mrproper: clean
	@mv  rapport.pdf $(DIR_TRASH)
