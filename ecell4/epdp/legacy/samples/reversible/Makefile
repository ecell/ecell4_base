LOGLEVEL = ERROR
PYTHONPATH = ../..
PYTHON = python

OUTFILES=rev.3.out rev.2.out rev.1.out rev.0.out rev.-1.out rev.-2.out #rev.-3.out

all: $(OUTFILES)

clean:
	rm -f $(OUTFILES)

rev.3.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-2 80000000

rev.2.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-3 40000000

rev.1.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-4 20000000

rev.0.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-5 10000000

rev.-1.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-6 10000000

rev.-2.out:
	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-7 10000000

#rev.-3.out:
#	env PYTHONPATH=$(PYTHONPATH) LOGLEVEL=$(LOGLEVEL) $(PYTHON) -O run.py $@ 1.25e-8 1000000
