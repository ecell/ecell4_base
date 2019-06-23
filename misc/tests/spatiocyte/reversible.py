# -*- coding: utf_8 -*-

import numpy as np
from ecell4 import *
from ecell4.extra.ensemble import ensemble_simulations

def main():
    radius, D = 5.0e-3, 1.0
    N_A = 60
    U = 0.5
    ka_factor = 0.1

    number_of_samples = 20

    kD = 4 * np.pi * (radius * 2) * (D * 2)
    ka = kD * ka_factor
    kd = ka * N_A * U * U / (1 - U)
    kon = ka * kD / (ka + kD)
    koff = kd * kon / ka

    with species_attributes():
        A | B | C | {'radius': str(radius), 'D': str(D)}

    with reaction_rules():
        A + B == C | (kon, koff)

    m = get_model()

    rng = GSLRandomNumberGenerator()
    rng.seed(0)

    y0 = {'A': N_A, 'B': N_A}
    duration = 3
    T = np.linspace(0, duration, 21)

    obs = run_simulation(np.linspace(0, duration, 101), y0,
                         model=ode.ODENetworkModel(m),
                         return_type='observer',
                         solver='ode')

    with species_attributes():
        A | B | C | {'radius': str(radius), 'D': str(D)}

    with reaction_rules():
        A + B == C | (ka, kd)

    m = get_model()

    ensemble_simulations(T, y0, model=m, return_type='matplotlib',
                         opt_args=('o', obs, '-'),
                         solver=('spatiocyte', radius),
                         n=number_of_samples)

if __name__ == '__main__':
    main()
