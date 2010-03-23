PYTHONPATH=../..

OUTFILES=\
out-C.py\
out-V.py\
out-N300.py\
out-N3000.py


all: $(OUTFILES)

out-C.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py C

out-V.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py V

out-N300.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py N300

out-N3000.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py N3000

clean:
	rm -f $(OUTFILES)