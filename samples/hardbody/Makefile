PYTHONPATH=../..

OUTFILES=\
out_C.py\
out_V.py\
out_N300.py\
out_N3000.py


all: $(OUTFILES)

out_C.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py C

out_V.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py V

out_N300.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py N300

out_N3000.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py N3000

clean:
	rm -f $(OUTFILES)