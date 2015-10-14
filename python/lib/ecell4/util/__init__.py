from .decorator import reaction_rules, species_attributes, parameters, get_model, reset_model
from . import viz

__all__ = [
    'run_simulation', 'ensemble_simulations', 'load_world',
    'reaction_rules', 'species_attributes', 'parameters', 'get_model', 'reset_model',
    'viz']


def load_world(filename):
    """
    Load a world from the given HDF5 filename.
    The return type is determined by ``ecell4.core.load_version_information``.

    Parameters
    ----------
    filename : str
        A HDF5 filename.

    Returns
    -------
    w : World
        Return one from ``BDWorld``, ``EGFRDWorld``, ``MesoscopicWorld``,
        ``ODEWorld``, ``GillespieWorld`` and ``LatticeWorld``.

    """
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
        t, y0={}, volume=1.0, model=None, solver='ode',
        factory=None, is_netfree=False, species_list=None, without_reset=False,
        return_type='matplotlib', plot_args=(), plot_kwargs={}):
    """Run a simulation with the given model and plot the result on IPython
    notebook with matplotlib.

    Parameters
    ----------
    t : array
        A sequence of time points for which to solve for 'm'.
    y0 : dict
        Initial condition.
    volume : Real, optional
    model : Model, optional
    solver : str, optional
        Solver type. Choose one from 'ode', 'gillespie', 'lattice', 'meso',
        'bd' and 'egfrd'. Default is 'ode'.
    species_list : list of str, optional
        A list of names of Species observed. If None, log all.
        Default is None.
    return_type : str, optional
        Choose a type of return value from 'array', 'observer',
        'matplotlib', 'nyaplot' or None.
        If None, return and plot nothing. Default is 'matplotlib'.
    plot_args: list, tuple or dict, optional
        Arguments for plotting. If return_type suggests no plotting, just ignored.
    plot_kwargs: dict, optional
        Arguments for plotting. If return_type suggests no plotting or
        plot_args is a list or tuple, just ignored.
        i.e.) viz.plot_number_observer(obs, *plot_args, **plot_kwargs)
    factory: Factory, optional
    is_netfree: bool, optional
        Whether the model is netfree or not. When a model is given as an
        argument, just ignored. Default is False.

    Returns
    -------
    value : list, TimingNumberObserver, or None
        Return a value suggested by ``return_type``.
        When ``return_type`` is 'array', return a time course data.
        When ``return_type`` is 'observer', return an observer.
        Return nothing if else.

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

    if isinstance(w, ecell4.ode.ODEWorld):
        # w.bind_to(model)  # stop binding for ode
        for serial, n in y0.items():
            w.set_value(ecell4.Species(serial), n)
    else:
        w.bind_to(model)
        for serial, n in y0.items():
            w.add_molecules(ecell4.Species(serial), n)

    if species_list is None:
        if isinstance(model, ecell4.ode.ODENetworkModel):
            #XXX: A bit messy way
            species_list = [sp.serial() for sp in model.list_species()]
        else:
            seeds = [ecell4.Species(serial) for serial in y0.keys()]
            species_list = [
                sp.serial() for sp in model.expand(seeds).list_species()]

    obs = ecell4.TimingNumberObserver(t, species_list)
    sim = f.create_simulator(model, w)
    # sim = f.create_simulator(w)
    sim.run(t[-1], obs)

    if return_type == 'matplotlib':
        if isinstance(plot_args, (list, tuple)):
            ecell4.viz.plot_number_observer(obs, *plot_args, **plot_kwargs)
        elif isinstance(plot_args, dict):
            # plot_kwargs is ignored
            ecell4.viz.plot_number_observer(obs, **plot_args)
        else:
            raise ValueError('plot_args [{}] must be list or dict.'.format(
                repr(plot_args)))
    elif return_type == 'nyaplot':
        if isinstance(plot_args, list):
            ecell4.viz.plot_number_observer_with_nya(obs, *plot_args, **plot_kwargs)
        elif isinstance(plot_args, dict):
            # plot_kwargs is ignored
            ecell4.viz.plot_number_observer_with_nya(obs, **plot_args)
        else:
            raise ValueError('plot_args [{}] must be list or dict.'.format(
                repr(plot_args)))
    elif return_type == 'observer':
        return obs
    elif return_type == 'array':
        return obs.data()

def ensemble_simulations(N=1, *args, **kwargs):
    """Run simulations with the given model and take the ensemble.

    Parameters
    ----------
    N : int
        A number of samples.
    args : list
    kwargs : dict
        Arguments for ``run_simulation``.
        See the help of ``run_simulation`` for details.

    Returns
    -------
    value : list, DummyObserver, or None
        Return a value suggested by ``return_type``.
        When ``return_type`` is 'array', return a time course data.
        When ``return_type`` is 'observer', return an observer.
        Return nothing if else.

    See Also
    --------
    run_simulation

    """
    import ecell4

    return_type = kwargs.get('return_type', 'matplotlib')
    plot_args = kwargs.get('plot_args', ())
    plot_kwargs = kwargs.get('plot_kwargs', {})
    kwargs.update({'return_type': 'observer'})

    import numpy

    class DummyObserver(object):

        def __init__(self, targets, data, func=lambda x: numpy.asarray(x)):
            self.__targets = targets
            self.__func = func
            self.__data = self.__func(data)

        def targets(self):
            return self.__targets

        def data(self):
            return self.__data.copy()

        def append(self, data):
            self.__data += self.__func(data)

    tmp = ecell4.util.run_simulation(*args, **kwargs)
    obs = DummyObserver(tmp.targets(), tmp.data(), lambda x: numpy.array(x, numpy.float64) / N)
    # import time
    for i in range(N - 1):
        # time.sleep(1)  #XXX: This is a missy way to vary random number seeds.
        tmp = ecell4.util.run_simulation(*args, **kwargs)
        obs.append(tmp.data())

    if return_type == 'matplotlib':
        if isinstance(plot_args, (list, tuple)):
            ecell4.viz.plot_number_observer(obs, *plot_args, **plot_kwargs)
        elif isinstance(plot_args, dict):
            # plot_kwargs is ignored
            ecell4.viz.plot_number_observer(obs, **plot_args)
        else:
            raise ValueError('plot_args [{}] must be list or dict.'.format(
                repr(plot_args)))
    elif return_type == 'nyaplot':
        if isinstance(plot_args, list):
            ecell4.viz.plot_number_observer_with_nya(obs, *plot_args, **plot_kwargs)
        elif isinstance(plot_args, dict):
            # plot_kwargs is ignored
            ecell4.viz.plot_number_observer_with_nya(obs, **plot_args)
        else:
            raise ValueError('plot_args [{}] must be list or dict.'.format(
                repr(plot_args)))
    elif return_type == 'observer':
        return obs
    elif return_type == 'array':
        return obs.data()
