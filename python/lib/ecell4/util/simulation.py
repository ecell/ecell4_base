import collections

from .decorator import get_model, reset_model
from . import viz


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
        ``ODEWorld``, ``GillespieWorld`` and ``SpatiocyteWorld``.

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
    elif vinfo.startswith("ecell4-spatiocyte"):
        return ecell4.spatiocyte.SpatiocyteWorld(filename)
    elif vinfo == "":
        raise RuntimeError("No version information was found in [{0}]".format(filename))
    raise RuntimeError("Unkown version information [{0}]".format(vinfo))

def get_factory(solver, *args):
    import ecell4

    if solver == 'ode':
        return ecell4.ode.ODEFactory(*args)
    elif solver == 'gillespie':
        return ecell4.gillespie.GillespieFactory(*args)
    elif solver == 'spatiocyte':
        return ecell4.spatiocyte.SpatiocyteFactory(*args)
    elif solver == 'meso':
        return ecell4.meso.MesoscopicFactory(*args)
    elif solver == 'bd':
        return ecell4.bd.BDFactory(*args)
    elif solver == 'egfrd':
        return ecell4.egfrd.EGFRDFactory(*args)
    else:
        raise ValueError(
            'unknown solver name was given: ' + repr(solver)
            + '. use ode, gillespie, spatiocyte, meso, bd or egfrd')

def list_species(model, seeds=None):
    seeds = None or []

    from ecell4.ode import ODENetworkModel
    from ecell4 import Species
    if isinstance(model, ODENetworkModel):
        #XXX: A bit messy way
        return sorted([sp.serial() for sp in model.list_species()])

    if not isinstance(seeds, list):
        seeds = list(seeds)

    expanded = model.expand([Species(serial) for serial in seeds])
    species_list = [sp.serial() for sp in expanded.list_species()]
    species_list = sorted(set(seeds + species_list))
    return species_list

def run_simulation(
        t, y0=None, volume=1.0, model=None, solver='ode',
        is_netfree=False, species_list=None, without_reset=False,
        return_type='matplotlib', opt_args=(), opt_kwargs=None,
        structures=None, observers=(), progressbar=0, rndseed=None,
        factory=None, ## deprecated
        **kwargs):
    """Run a simulation with the given model and plot the result on IPython
    notebook with matplotlib.

    Parameters
    ----------
    t : array or Real
        A sequence of time points for which to solve for 'm'.
    y0 : dict
        Initial condition.
    volume : Real or Real3, optional
        A size of the simulation volume.
        Keyword 'v' is a shortcut for specifying 'volume'.
    model : Model, optional
        Keyword 'm' is a shortcut for specifying 'model'.
    solver : str, tuple or Factory, optional
        Solver type. Choose one from 'ode', 'gillespie', 'spatiocyte', 'meso',
        'bd' and 'egfrd'. Default is 'ode'.
        When tuple is given, the first value must be str as explained above.
        All the rest is used as arguments for the corresponding factory class.
        Keyword 's' is a shortcut for specifying 'solver'.
    species_list : list of str, optional
        A list of names of Species observed. If None, log all.
        Default is None.
    return_type : str, optional
        Choose a type of return value from 'array', 'observer',
        'matplotlib', 'nyaplot', 'world', 'dataframe' or None.
        If None, return and plot nothing. Default is 'matplotlib'.
        'dataframe' requires numpy and pandas libraries.
        Keyword 'r' is a shortcut for specifying 'return_type'.
    opt_args: list, tuple or dict, optional
        Arguments for plotting. If return_type suggests no plotting, just ignored.
    opt_kwargs: dict, optional
        Arguments for plotting. If return_type suggests no plotting or
        opt_args is a list or tuple, just ignored.
        i.e.) viz.plot_number_observer(obs, *opt_args, **opt_kwargs)
    is_netfree: bool, optional
        Whether the model is netfree or not. When a model is given as an
        argument, just ignored. Default is False.
    structures : dict, optional
        A dictionary which gives pairs of a name and shape of structures.
        Not fully supported yet.
    observers : Observer or list, optional
        A list of extra observer references.
    progressbar : float, optional
        A timeout for a progress bar in seconds.
        When the value is not more than 0, show nothing.
        Default is 0.
    rndseed : int, optional
        A random seed for a simulation.
        This argument will be ignored when 'solver' is given NOT as a string.

    Returns
    -------
    value : list, TimingNumberObserver, World or None
        Return a value suggested by ``return_type``.
        When ``return_type`` is 'array', return a time course data.
        When ``return_type`` is 'observer', return an observer.
        When ``return_type`` is 'world', return the last state of ``World``.
        Return nothing if else.

    """
    y0 = y0 or {}
    opt_kwargs = opt_kwargs or {}
    structures = structures or {}

    for key, value in kwargs.items():
        if key == 'r':
            return_type = value
        elif key == 'v':
            volume = value
        elif key == 's':
            solver = value
        elif key == 'm':
            model = value
        else:
            raise ValueError(
                "An unknown keyword argument was given [{}={}]".format(key, value))

    import ecell4

    if factory is not None:
        # f = factory  #XXX: will be deprecated in the future. just use solver
        raise ValueError(
            "Argument 'factory' is no longer available. Use 'solver' instead.")
    elif isinstance(solver, str):
        f = get_factory(solver)
    elif isinstance(solver, collections.Iterable):
        f = get_factory(*solver)
    else:
        f = solver

    if rndseed is not None:
        f = f.rng(ecell4.GSLRandomNumberGenerator(rndseed))

    if model is None:
        model = ecell4.util.decorator.get_model(is_netfree, without_reset)

    if isinstance(volume, ecell4.Real3):
        edge_lengths = volume
    else:
        L = ecell4.cbrt(volume)
        edge_lengths = ecell4.Real3(L, L, L)

    w = f.create_world(edge_lengths)

    for (name, shape) in structures.items():
        w.add_structure(ecell4.Species(name), shape)

    if isinstance(w, ecell4.ode.ODEWorld):
        # w.bind_to(model)  # stop binding for ode
        for serial, n in y0.items():
            w.set_value(ecell4.Species(serial), n)
    else:
        w.bind_to(model)
        for serial, n in y0.items():
            w.add_molecules(ecell4.Species(serial), n)

    if species_list is None:
        species_list = list_species(model, y0.keys())

    if not isinstance(t, collections.Iterable):
        t = [float(t) * i / 100 for i in range(101)]

    obs = ecell4.TimingNumberObserver(t, species_list)
    sim = f.create_simulator(model, w)
    # sim = f.create_simulator(w)

    if not isinstance(observers, collections.Iterable):
        observers = (observers, )
    if return_type not in ('world', None):
        observers = (obs, ) + tuple(observers)

    if progressbar > 0:
        from .progressbar import progressbar as pb
        pb(sim, timeout=progressbar, flush=True).run(t[-1], observers)
    else:
        sim.run(t[-1], observers)

    if return_type in ('matplotlib', 'm'):
        if isinstance(opt_args, (list, tuple)):
            ecell4.viz.plot_number_observer(obs, *opt_args, **opt_kwargs)
        elif isinstance(opt_args, dict):
            # opt_kwargs is ignored
            ecell4.viz.plot_number_observer(obs, **opt_args)
        else:
            raise ValueError('opt_args [{}] must be list or dict.'.format(
                repr(opt_args)))
    elif return_type in ('nyaplot', 'n'):
        if isinstance(opt_args, (list, tuple)):
            ecell4.viz.plot_number_observer_with_nya(obs, *opt_args, **opt_kwargs)
        elif isinstance(opt_args, dict):
            # opt_kwargs is ignored
            ecell4.viz.plot_number_observer_with_nya(obs, **opt_args)
        else:
            raise ValueError('opt_args [{}] must be list or dict.'.format(
                repr(opt_args)))
    elif return_type in ('observer', 'o'):
        return obs
    elif return_type in ('array', 'a'):
        return obs.data()
    elif return_type in ('dataframe', 'd'):
        import pandas
        import numpy
        data = numpy.array(obs.data()).T
        return pandas.concat([
            pandas.DataFrame(dict(Time=data[0], Value=data[i + 1],
                                  Species=sp.serial(), **opt_kwargs))
            for i, sp in enumerate(obs.targets())])
    elif return_type in ('world', 'w'):
        return sim.world()

def ensemble_simulations(N=1, *args, **kwargs):
    """Deprecated"""
    raise RuntimeError(
        "This function was deprecated. Use ecell4.extra.ensemble.ensemble_simulations instread.")
