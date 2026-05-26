# Makefile minimo para compilar las versiones inglesa y espanola con latexmk.

MAIN = main
MAIN_ES = main_es
OUTDIR = build/pdf

.PHONY: build build-en build-es watch watch-en watch-es clean clean-en clean-es

build: build-en

build-en:
	@mkdir -p "$(OUTDIR)"
	latexmk -pdf -outdir="$(OUTDIR)" "$(MAIN).tex"

build-es:
	@mkdir -p "$(OUTDIR)"
	latexmk -pdf -outdir="$(OUTDIR)" "$(MAIN_ES).tex"

watch: watch-en

watch-en:
	@mkdir -p "$(OUTDIR)"
	latexmk -pvc -pdf -outdir="$(OUTDIR)" "$(MAIN).tex"

watch-es:
	@mkdir -p "$(OUTDIR)"
	latexmk -pvc -pdf -outdir="$(OUTDIR)" "$(MAIN_ES).tex"

clean: clean-en clean-es

clean-en:
	latexmk -c -outdir="$(OUTDIR)" "$(MAIN).tex"

clean-es:
	latexmk -c -outdir="$(OUTDIR)" "$(MAIN_ES).tex"
