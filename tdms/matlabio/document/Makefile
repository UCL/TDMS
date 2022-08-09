.PHONY: pdf clean

BUILD_DIR = build
DOC_FILE = matlab_array.pdf

pdf: $(DOC_FILE)

$(DOC_FILE): matlab_array.tex
	mkdir -p $(BUILD_DIR)
	pdflatex -synctex=1 -interaction=nonstopmode -output-directory=$(BUILD_DIR) $<
	cp $(BUILD_DIR)/$@ $@

clean:
	rm -rf $(BUILD_DIR) $(DOC_FILE)
