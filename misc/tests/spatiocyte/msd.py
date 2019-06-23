# -*- coding: utf_8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
from ecell4 import *

def main():
    radius, D = 5.0e-3, 1.0
    model = NetworkModel()
    model.add_species_attribute(Species("A", str(radius), str(D)))

    rng = GSLRandomNumberGenerator()
    rng.seed(0)

    factory = spatiocyte.SpatiocyteFactory(radius).rng(rng) 
    def calc_squared_displacements(trajectory):
        origin = trajectory[0]
        return list(map(lambda pos: length_sq(pos - origin), trajectory))

    def run_and_calc_msd(duration):
        world = factory.create_world(Real3(1.0, 1.0, 1.0))
        world.bind_to(model)
        world.add_molecules(Species("A"), 60)

        obs = FixedIntervalTrajectoryObserver(0.01)
        simulator = factory.create_simulator(world)
        simulator.run(duration, obs)

        times = np.array(obs.t())
        msds = np.mean(list(map(calc_squared_displacements, obs.data())), axis=0)

        return times, msds

    def test_msd(num, duration):
        times, msds = run_and_calc_msd(duration)
        for _ in range(0, num - 1):
            msds += run_and_calc_msd(duration)[1]
        msds /= num
        return times, msds

    times, msds = test_msd(10, 1.00)


    # Plot

    plt.plot(times, 6*D*times, 'k-', label='Expected')
    plt.plot(times[::10], msds[::10], 'ro', label='Spatiocyte')
    plt.xlabel('Time')
    plt.ylabel('Mean Squared Displacement')
    plt.legend(loc='best')
    plt.title('MSD in a space')

    name=os.path.splitext(os.path.basename(__file__))[0]
    plt.savefig(name+'.png')

if __name__ == '__main__':
    main()
