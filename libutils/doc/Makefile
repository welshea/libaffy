
####
## 
## Makefile for doc
##
## (invoked from the top-level libutils makefile)
##
####

## All user-defineable parameters are set at the top-level; you should
## not need to modify this makefile.

MAKEINFO?=makeinfo
TEXI2PDF?=texi2pdf
RM?=rm

JUNK_EXTS:= aux cp fn fns ky log pg toc tp vr vrs

BASENAME:= libutils

SOURCE:= $(BASENAME).texi 

.PHONY: all
all: html pdf

.PHONY: pdf
pdf: $(BASENAME).pdf

.PHONY: html
html:   $(SOURCE)
	$(MAKEINFO) --html -ohtml $(SOURCE)
	$(MAKEINFO) --html -o$(BASENAME).html --no-split $(SOURCE)

$(BASENAME).pdf: $(SOURCE)
	$(TEXI2PDF) $(SOURCE)

.PHONY: clean
clean:
	-$(RM) -fr html $(BASENAME).pdf $(BASENAME).html
	-$(RM) $(foreach x,$(JUNK_EXTS),*.$(x))
