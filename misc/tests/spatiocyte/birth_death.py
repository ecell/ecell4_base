# -*- coding: utf_8 -*-

import numpy as np
from ecell4 import *
from ecell4.extra.ensemble import ensemble_simulations

radius, D = 5.0e-3, 1.0

with species_attributes():
    A | {'radius': str(radius), 'D': str(D)}

with reaction_rules():
    ~A > A | 45.0
    A > ~A | 1.5

model = get_model()

rng = GSLRandomNumberGenerator()
rng.seed(0)

number_of_samples = 20
y0 = {}
duration = 3
T = np.linspace(0, duration, 21)
V = 8

obs = run_simulation(np.linspace(0, duration, 101), y0, volume=V,
                     model=ode.ODENetworkModel(model),
                     return_type='observer',
                     solver='ode')

ensemble_simulations(T, y0, volume=V, model=model,
                     return_type='matplotlib',
                     opt_args=('o', obs, '-'),
                     solver=('spatiocyte', radius),
                     n=number_of_samples)
