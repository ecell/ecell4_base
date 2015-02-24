from decorator import reaction_rules, species_attributes, get_model, reset_model
import viz

__all__ = [
    'run_simulation', 'reaction_rules', 'species_attributes',
    'get_model', 'reset_model',
    'viz']


def run_simulation(
        t, y0, volume=1.0, model=None, with_plot=True, solver='ode',
        is_netfree=False, species_list=None):
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
    is_netfree: bool, optional
        Whether the model is netfree or not. When a model is given as an
        argument, just ignored.
    """
    import ecell4

    if solver == 'ode':
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
        model = ecell4.util.decorator.get_model(is_netfree, True)

    edge_lengths = ecell4.Real3(L, L, L)
    w = f.create_world(edge_lengths)
    w.bind_to(model)
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
        ecell4.viz.plot_number_observer(obs)

    return obs.data()
