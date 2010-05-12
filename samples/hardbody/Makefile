PYTHONPATH=../..

OUTFILES=\
out_C.py\
out_V.py\
out_N300.py\
out_N3000.py\
out_BD.py\
out_BD2.py


all: $(OUTFILES)

out_C.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py C

out_V.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py V

out_N300.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py N300

out_N3000.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py N3000

out_BD.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py BD

out_BD2.py:
	PYTHONPATH=$(PYTHONPATH) python -O run_all.py BD2

clean:
	rm -f $(OUTFILES)