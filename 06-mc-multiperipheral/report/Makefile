# Build mumu_cancellation.pdf (figures first, then two pdflatex passes).
# PDF is local-only on main; publish via hooks/publish-pdf-builds.sh.

TEX       = mumu_cancellation.tex
PDF       = mumu_cancellation.pdf
FIG_TEX   = $(wildcard figures/*.tex)
FIG_PDF   = $(FIG_TEX:.tex=.pdf)
LATEX     = pdflatex
LATEXFLAGS = -interaction=nonstopmode -halt-on-error

.PHONY: all pdf figures clean distclean pdf-builds

all: pdf

pdf: $(PDF)

figures: $(FIG_PDF)

figures/%.pdf: figures/%.tex
	cd figures && $(LATEX) $(LATEXFLAGS) $(notdir $<)

$(PDF): $(TEX) $(FIG_PDF)
	$(LATEX) $(LATEXFLAGS) $(TEX)
	$(LATEX) $(LATEXFLAGS) $(TEX)

clean:
	rm -f *.aux *.log *.out *.toc *.fls *.fdb_latexmk *.synctex.gz
	rm -f figures/*.aux figures/*.log

distclean: clean
	rm -f $(PDF) $(FIG_PDF)

# One-shot: build and force-push orphan pdf-builds (both remotes).
pdf-builds:
	../../hooks/publish-pdf-builds.sh
