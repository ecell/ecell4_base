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

def get_factory(solver, *args):
    import ecell4

    if solver == 'ode':
        return ecell4.ode.ODEFactory(*args)
    elif solver == 'gillespie':
        return ecell4.gillespie.GillespieFactory(*args)
    elif solver == 'lattice':
        return ecell4.lattice.LatticeFactory(*args)
    elif solver == 'meso':
        return ecell4.meso.MesoscopicFactory(*args)
    elif solver == 'bd':
        return ecell4.bd.BDFactory(*args)
    elif solver == 'egfrd':
        return ecell4.egfrd.EGFRDFactory(*args)
    else:
        raise ValueError(
            'unknown solver name was given: ' + repr(solver)
            + '. use ode, gillespie, lattice, meso, bd or egfrd')

def run_simulation(
        t, y0={}, volume=1.0, model=None, solver='ode',
        factory=None, is_netfree=False, species_list=None, without_reset=False,
        return_type='matplotlib', plot_args=(), plot_kwargs={},
        structures={}, observers=()):
    """Run a simulation with the given model and plot the result on IPython
    notebook with matplotlib.

    Parameters
    ----------
    t : array
        A sequence of time points for which to solve for 'm'.
    y0 : dict
        Initial condition.
    volume : Real or Real3, optional
        A size of the simulation volume.
    model : Model, optional
    solver : str, optional
        Solver type. Choose one from 'ode', 'gillespie', 'lattice', 'meso',
        'bd' and 'egfrd'. Default is 'ode'.
    species_list : list of str, optional
        A list of names of Species observed. If None, log all.
        Default is None.
    return_type : str, optional
        Choose a type of return value from 'array', 'observer',
        'matplotlib', 'nyaplot', 'world' or None.
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
    structures : dict, optional
        A dictionary which gives pairs of a name and shape of structures.
        Not fully supported yet.
    observers : list, optional
        A list of extra observer references.

    Returns
    -------
    value : list, TimingNumberObserver, World or None
        Return a value suggested by ``return_type``.
        When ``return_type`` is 'array', return a time course data.
        When ``return_type`` is 'observer', return an observer.
        When ``return_type`` is 'world', return the last state of ``World``.
        Return nothing if else.

    """
    import ecell4

    if factory is not None:
        f = factory
    else:
        f = get_factory(solver)

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
    sim.run(t[-1], (obs, ) + tuple(observers))

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
    elif return_type == 'world':
        return sim.world()

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

    # if (len(args) <= 5 or args[5] is None) and 'factory' not in kwargs.keys():
    #     if len(args) <= 4 and 'solver' not in kwargs.keys():
    #         kwargs['factory'] = ecell4.ode.ODEFactory()
    #     else:
    #         solver = kwargs['solver'] if len(args) <= 4 else args[4]
    #         if solver == 'ode':
    #             kwargs['factory'] = get_factory(solver)
    #         else:
    #             kwargs['factory'] = get_factory(
    #                 solver, ecell4.core.GSLRandomNumberGenerator())
    #         # 'solver' would be ignored in run_simulation

    errorbar = kwargs.pop('errorbar', True)

    return_type = kwargs.get('return_type', 'matplotlib') if len(args) <= 9 else args[9]
    if return_type == 'world':
        raise ValueError('return_type "world" is not supported in ensemble_simulations.')
    if len(args) > 9:
        args[9] = 'observer'
    else:
        kwargs.update({'return_type': 'observer'})

    plot_args = kwargs.get('plot_args', ()) if len(args) <= 10 else args[10]
    plot_kwargs = kwargs.get('plot_kwargs', {}) if len(args) <= 11 else args[11]

    import numpy

    class DummyObserver(object):

        def __init__(self, targets, data):
            self.__targets = targets
            self.__t = numpy.array(data).T[0]
            self.__tot = self.trancate(data)
            self.__counts = 1

        def targets(self):
            return self.__targets

        def t(self):
            return self.__t

        def counts(self):
            return self.__counts

        def trancate(self, data):
            return numpy.array(data, numpy.float64).T[1: ]

        def average_(self):
            return self.__tot / self.__counts

        def average(self):
            return numpy.vstack([self.__t, self.average_()]).T

        def data(self):
            return self.average()

        def append(self, data):
            self.__tot += self.trancate(data)
            self.__counts += 1

    class DummyObserverWithError(DummyObserver):

        def __init__(self, targets, data):
            DummyObserver.__init__(self, targets, data)
            self.__sqtot = self.trancate(data) ** 2

        def std_(self):
            return self.__sqtot / self.counts() - self.average_() ** 2

        def error_(self):
            return self.std_() / numpy.sqrt(self.counts())

        def error(self):
            return numpy.vstack([self.t(), self.error_()]).T

        def append(self, data):
            DummyObserver.append(self, data)
            self.__sqtot += self.trancate(data) ** 2

    tmp = ecell4.util.run_simulation(*args, **kwargs)

    if errorbar:
        obs = DummyObserverWithError(tmp.targets(), tmp.data())
    else:
        obs = DummyObserver(tmp.targets(), tmp.data())

    for i in range(N - 1):
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
