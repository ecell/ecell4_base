
# BUILT
if HAVE_PDFLATEX
# Automake manual 18.2 Handling new file extensions
.tex.pdf:
	pdflatex $<

# Compile and install; do not distribute.
pdf_DATA = \
	implementation_notes.pdf\
	p1_fp.pdf
endif


# INSTALL
# Install and distribute.
nobase_dist_doc_DATA = \
	testing.txt\
	math/*


# DIST
# Distribute; do not install.
EXTRA_DIST = *.tex


# CLEAN
clean-local:
	rm -f *.log *.aux *.pdf

