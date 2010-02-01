LOGLEVEL = ERROR
PYTHONPATH = ../..
PYTHON = python

all: rev.3.out rev.2.out rev.1.out rev.0.out rev.-1.out rev.-2.out rev.-3.out

rev.3.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py rev.3.out 1.25e-2 10000000

rev.2.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py rev.2.out 1.25e-3 10000000

rev.1.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py rev.1.out 1.25e-4 4000000

rev.0.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py rev.0.out 1.25e-5 4000000

rev.-1.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py rev.-1.out 1.25e-6 4000000

rev.-2.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py rev.-2.out 1.25e-7 4000000

rev.-3.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py rev.-3.out 1.25e-8 1000000
