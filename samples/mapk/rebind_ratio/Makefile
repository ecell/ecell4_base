



DATADIR=/home/shafi/run/mapk/09/data/

PATTERNS=\
mapk3_1e-15_0.25_fixed_0_normal_ALL_reactions.dat\
mapk3_1e-15_0.5_fixed_0_normal_ALL_reactions.dat\
mapk3_1e-15_1_fixed_0_normal_ALL_reactions.dat\
mapk3_1e-15_2_fixed_0_normal_ALL_reactions.dat\
mapk3_1e-15_4_fixed_0_normal_ALL_reactions.dat\
mapk3_1e-15_0.25_fixed_1e-6_normal_ALL_reactions.dat\
mapk3_1e-15_0.5_fixed_1e-6_normal_ALL_reactions.dat\
mapk3_1e-15_1_fixed_1e-6_normal_ALL_reactions.dat\
mapk3_1e-15_2_fixed_1e-6_normal_ALL_reactions.dat\
mapk3_1e-15_4_fixed_1e-6_normal_ALL_reactions.dat\
mapk3_1e-15_0.25_fixed_1e-4_normal_ALL_reactions.dat\
mapk3_1e-15_0.5_fixed_1e-4_normal_ALL_reactions.dat\
mapk3_1e-15_1_fixed_1e-4_normal_ALL_reactions.dat\
mapk3_1e-15_2_fixed_1e-4_normal_ALL_reactions.dat\
mapk3_1e-15_4_fixed_1e-4_normal_ALL_reactions.dat\
mapk3_1e-15_0.25_fixed_1e-2_normal_ALL_reactions.dat\
mapk3_1e-15_0.5_fixed_1e-2_normal_ALL_reactions.dat\
mapk3_1e-15_1_fixed_1e-2_normal_ALL_reactions.dat\
mapk3_1e-15_2_fixed_1e-2_normal_ALL_reactions.dat\
mapk3_1e-15_4_fixed_1e-2_normal_ALL_reactions.dat\
mapk3_1e-15_0.25_fixed_1e-1_normal_ALL_reactions.dat\
mapk3_1e-15_0.5_fixed_1e-1_normal_ALL_reactions.dat\
mapk3_1e-15_1_fixed_1e-1_normal_ALL_reactions.dat\
mapk3_1e-15_2_fixed_1e-1_normal_ALL_reactions.dat\
mapk3_1e-15_4_fixed_1e-1_normal_ALL_reactions.dat\

REBIND=$(PATTERNS:.dat=.rebind)


$(REBIND):
	python ../rebind_ratio.py $(DATADIR)/$(@:.rebind=.dat); \


ratio: $(REBIND)
	for i in $(REBIND); do\
		echo $$i; \
		python ../../../bin/event_ratio.py $$i; \
	done

plot:
	python plot.py



clean:
	rm -f *.rebind
