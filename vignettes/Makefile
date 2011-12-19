All: crlmmDownstream clean

# The first target must create all the pdfs
crlmmDownstream: crlmmDownstream.tex
	cp ../scripts/crlmmDownstream.pdf .
	texi2dvi --pdf VanillaICE.tex

clean:
	-$(RM) -f *.out *.bbl *.log *.aux *.blg *.brf *.toc *.tex