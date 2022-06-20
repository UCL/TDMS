.PHONY: pdf clean

BUILD_DIR = build
DOC_FILE = masterdoc.pdf

pdf: $(DOC_FILE)

$(DOC_FILE): masterdoc.tex refs.bib $(wildcard figures/*.pdf)
	mkdir -p $(BUILD_DIR)
	pdflatex -synctex=1 -interaction=nonstopmode -output-directory=$(BUILD_DIR) $<
	bibtex $(BUILD_DIR)/$(basename $@)
	pdflatex -synctex=1 -interaction=nonstopmode -output-directory=$(BUILD_DIR) $<
	pdflatex -synctex=1 -interaction=nonstopmode -output-directory=$(BUILD_DIR) $<
	cp $(BUILD_DIR)/$@ $@

clean:
	rm -rf $(BUILD_DIR) $(DOC_FILE)
