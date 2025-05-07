# From the latexmk Users' Guide (or 'man latexmk')

export TEXINPUTS=:../Common_Guide_Files//:
export BIBINPUTS=:../Common_Guide_Files/

.PHONY : FORCE_MAKE
TARGET:=Xyce_RG

all : $(TARGET).pdf
%.pdf : %.tex FORCE_MAKE
	latexmk -pdf $<

clean:
	latexmk -c -bibtex

realclean:
	latexmk -C -pdf -bibtex

outside: $(TARGET)_SANDoutside.pdf
$(TARGET)_SANDoutside.pdf: $(TARGET)_SANDoutside.tex
	latexmk -f -pdf $(TARGET)_SANDoutside.tex


remove-style: $(TARGET)_SANDoutside.tex
$(TARGET)_SANDoutside.tex: $(TARGET).tex
	cp $(TARGET).tex $(TARGET)_SANDoutside.tex && \
		sed -i \
		's/\\documentclass\[11pt,report\]{SANDreport}/\\documentclass\[11pt,letterpaper\]{scrreprt}/g' \
		$(TARGET)_SANDoutside.tex && \
		sed -i \
		's/\\usepackage\[sand\]{optional}/\\usepackage\[report\]{optional}/g' \
		$(TARGET)_SANDoutside.tex && \
		sed -i 's/SANDauthor/author/g' \
		$(TARGET)_SANDoutside.tex

