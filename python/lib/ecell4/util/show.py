from __future__ import print_function

import ecell4

from .viz import plot_number_observer, plot_trajectory, plot_world
from .simulation import load_world


def dump_model(m):
    if isinstance(m, ecell4.ode.ODENetworkModel):
        res = []
    else:
        res = [sp.serial() + "|" + str(dict(sp.list_attributes())) for sp in m.species_attributes()]
    res += [rr.as_string() for rr in m.reaction_rules()]
    # return res
    print('\n'.join(res))

def show(target, *args, **kwargs):
    """
    An utility function to display the given target object in the proper way.

    Paramters
    ---------
    target : NumberObserver, TrajectoryObserver, World, str
        When a NumberObserver object is given, show it with viz.plot_number_observer.
        When a TrajectoryObserver object is given, show it with viz.plot_trajectory_observer.
        When a World or a filename suggesting HDF5 is given, show it with viz.plot_world.

    """
    if isinstance(target, (ecell4.FixedIntervalNumberObserver, ecell4.NumberObserver, ecell4.TimingNumberObserver, )):
        plot_number_observer(target, *args, **kwargs)
    elif isinstance(target, (ecell4.FixedIntervalTrajectoryObserver, ecell4.FixedIntervalTrackingObserver)):
        plot_trajectory(target, *args, **kwargs)
    elif isinstance(target, (ecell4.ode.ODEWorld, ecell4.gillespie.GillespieWorld, ecell4.spatiocyte.SpatiocyteWorld, ecell4.meso.MesoscopicWorld, ecell4.bd.BDWorld, ecell4.egfrd.EGFRDWorld)):
        plot_world(target, *args, **kwargs)
    elif isinstance(target, (ecell4.core.Model, ecell4.core.NetworkModel, ecell4.core.NetfreeModel, ecell4.ode.ODENetworkModel)):
        dump_model(target)
    elif isinstance(target, str):
        try:
            w = simulation.load_world(target)
        except RuntimeError as e:
            raise ValueError("The given target [{}] is not supported.".format(repr(target)))
        else:
            show(w, *args, **kwargs)
    else:
        raise ValueError("The given target [{}] is not supported.".format(repr(target)))
