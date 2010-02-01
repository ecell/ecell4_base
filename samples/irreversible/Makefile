LOGLEVEL = ERROR
PYTHONPATH = ../..
PYTHON = python


all: irr.3.out irr.2.out irr.1.out irr.0.out irr.-1.out irr.-2.out irr.-3.out


irr.3.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py irr.3.out 1.25e-2 20000000

irr.2.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py irr.2.out 1.25e-3 20000000

irr.1.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py irr.1.out 1.25e-4 7000000

irr.0.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py irr.0.out 1.25e-5 5000000

irr.-1.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py irr.-1.out 1.25e-6 2000000

irr.-2.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py irr.-2.out 1.25e-7 2000000

irr.-3.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py irr.-3.out 1.25e-8 1000000
