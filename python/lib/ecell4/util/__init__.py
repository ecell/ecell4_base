from .decorator import reaction_rules, species_attributes, parameters, get_model, reset_model
from . import viz

__all__ = [
    'run_simulation', 'load_world',
    'reaction_rules', 'species_attributes', 'parameters', 'get_model', 'reset_model',
    'viz']


def load_world(filename):
    import ecell4

    vinfo = ecell4.core.load_version_information(filename)
    if vinfo.startswith("ecell4-bd"):
        return ecell4.bd.BDWorld(filename)
    elif vinfo.startswith("ecell4-egfrd"):
        return ecell4.egfrd.EGFRDWorld(filename)
    elif vinfo.startswith("ecell4-meso"):
        return ecell4.meso.MesoscopicWorld(filename)
    elif vinfo.startswith("ecell4-ode"):
        return ecell4.ode.ODEWorld(filename)
    elif vinfo.startswith("ecell4-gillespie"):
        return ecell4.gillespie.GillespieWorld(filename)
    elif vinfo.startswith("ecell4-lattice"):
        return ecell4.lattice.LatticeWorld(filename)
    elif vinfo == "":
        raise RuntimeError("No version information was found in [{0}]".format(filename))
    raise RuntimeError("Unkown version information [{0}]".format(vinfo))

def run_simulation(
        t, y0={}, volume=1.0, model=None, with_plot=True, solver='ode',
        factory=None, is_netfree=False, species_list=None, as_observer=False,
        without_reset=False):
    """Run a simulation with the given model and plot the result on IPython
    notebook with matplotlib.

    Parameters
    ----------
    t: array
        A sequence of time points for which to solve for 'm'.
    y0: dict
        Initial condition.
    volume: Real, optional
    model: Model, optional
    with_plot: bool, optional
        Whether to show the result as a plot.
    solver: str, optional
        Solver type. Choose one from 'ode', 'gillespie', 'lattice', 'meso',
        'bd' and 'egfrd'.
    factory: Factory, optional
    is_netfree: bool, optional
        Whether the model is netfree or not. When a model is given as an
        argument, just ignored.
    as_observer: bool, optional
        Return an Observer, but not an array.
    """
    import ecell4

    if factory is not None:
        f = factory
    elif solver == 'ode':
        f = ecell4.ode.ODEFactory()
    elif solver == 'gillespie':
        f = ecell4.gillespie.GillespieFactory()
    elif solver == 'lattice':
        f = ecell4.lattice.LatticeFactory()
    elif solver == 'meso':
        f = ecell4.meso.MesoscopicFactory()
    elif solver == 'bd':
        f = ecell4.bd.BDFactory()
    elif solver == 'egfrd':
        f = ecell4.egfrd.EGFRDFactory()
    else:
        raise ValueError(
            'unknown solver name was given: ' + repr(solver)
            + '. use ode, gillespie, lattice, meso, bd or egfrd')

    L = ecell4.cbrt(volume)

    if model is None:
        model = ecell4.util.decorator.get_model(is_netfree, without_reset)

    edge_lengths = ecell4.Real3(L, L, L)
    w = f.create_world(edge_lengths)
    w.bind_to(model)
    if isinstance(w, ecell4.ode.ODEWorld):
        for serial, n in y0.items():
            w.set_value(ecell4.Species(serial), n)
    else:
        for serial, n in y0.items():
            w.add_molecules(ecell4.Species(serial), n)

    if species_list is None:
        seeds = [ecell4.Species(serial) for serial in y0.keys()]
        species_list = [
            sp.serial() for sp in model.expand(seeds).list_species()]

    obs = ecell4.TimingNumberObserver(t, species_list)
    sim = f.create_simulator(w)
    sim.run(t[-1], obs)

    if with_plot:
        ecell4.viz.plot_number_observer_with_nya(obs)
        # ecell4.viz.plot_number_observer(obs)

    if as_observer:
        return obs

    return obs.data()
