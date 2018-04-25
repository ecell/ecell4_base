from ecell4.core import _save_bd5

def save_bd5(
        space, filename,
        group_index=0, object_name="molecule", spatial_unit="meter", time_unit="second",
        trunc=False, with_radius=False):
    """Save a space in the BDML-BD5 format (https://github.com/openssbd/BDML-BD5).

    Open file for read/write, if it already exists, and create a new file, otherwise.
    If trunc is True, always create a new file.
    A new group named `group_name` is created. If the group already exists, returns
    an exception.

    Parameters
    ----------
    space : Space
        A Space or World to be saved.
    filename : str
        A HDF5 filename.
    group_index : int, optional
        An index of the group written (0, 1, ..., n). Defaults to 0.
    object_name : str, optional
        A name of the object. Defaults to "molecule".
    spatial_unit : str, optional
        An unit of the length scale. Defaults to "meter".
    time_unit : str, optional
        An unit of the time scale. Defaults to "second".
    trunc : bool, optional
        Whether truncate file or not. If True, always overwrite the file when open it.
        Defaults to False.
    with_radius : bool, optional
        Defaults to False.

    Returns
    -------
    w : World
        Return one from ``BDWorld``, ``EGFRDWorld``, ``MesoscopicWorld``,
        ``ODEWorld``, ``GillespieWorld`` and ``SpatiocyteWorld``.
    """
    _save_bd5(space.as_base(), filename, group_index, object_name, spatial_unit, time_unit, trunc, with_radius)
