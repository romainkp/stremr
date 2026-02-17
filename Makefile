# Variables
INPUT = README.Rmd
HTML_OUT = README.html
MD_OUT = README.md

# Default target
all: $(HTML_OUT) $(MD_OUT)

# Rule to render both formats
$(HTML_OUT) $(MD_OUT): $(INPUT) inst/Bib/README-refs.bib
	Rscript -e "rmarkdown::render('$(INPUT)', output_format = 'all')"

# Clean up
clean:
	rm -f $(HTML_OUT) $(MD_OUT)
	rm -rf README_files
