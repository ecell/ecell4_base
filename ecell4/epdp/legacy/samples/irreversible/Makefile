LOGLEVEL = ERROR
PYTHONPATH = ../..
PYTHON = python

OUTFILES=irr.3.out irr.2.out irr.1.out irr.0.out irr.-1.out irr.-2.out #irr.-3.out

all: $(OUTFILES)

clean:
	rm -f $(OUTFILES)

irr.3.out:			
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-2 100000000

irr.2.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-3 70000000

irr.1.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-4 7000000

irr.0.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-5 5000000

irr.-1.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-6 5000000

irr.-2.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-7 5000000

#irr.-3.out:
#	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-8 1000000
