import pathlib

from ecell4.core import _save_bd5
from ecell4.util.simulation import load_world

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
    space : Space, str, pathlib.PurePath, list, tuple or set
        A Space or World to be saved. If str or pathlib.PurePath is given, a space is
        loaded from the given path. If this is an iterable (list, tuple, set), apply
        this function to each element of the given.
    filename : str
        A HDF5 filename.
    group_index : int, optional
        An index of the group written (0, 1, ..., n). Defaults to 0.
    object_name : str, optional
        A name of the object. Its length must be less than 128. Defaults to "molecule".
    spatial_unit : str, optional
        An unit of the length scale. Its length must be less than 16. Defaults to "meter".
    time_unit : str, optional
        An unit of the time scale. Its length must be less than 16. Defaults to "second".
    trunc : bool, optional
        Whether truncate file or not. If True, always overwrite the file when open it.
        Defaults to False.
    with_radius : bool, optional
        Whether save the radius of particles. If True, particles are saved as 'sphere',
        otherwise, as 'point'. Defaults to False.

    """
    if isinstance(space, (list, tuple, set)):
        for i, space_ in enumerate(space):
            assert not isinstance(space_, (list, tuple, set))
            save_bd5(
                space_, filename, group_index + i, object_name, spatial_unit, time_unit,
                trunc if i == 0 else False, with_radius)
    elif isinstance(space, str):
        save_bd5(load_world(space), filename, group_index, object_name, spatial_unit, time_unit, trunc, with_radius)
    elif isinstance(space, pathlib.PurePath):
        save_bd5(str(space), filename, group_index, object_name, spatial_unit, time_unit, trunc, with_radius)
    else:
        # space is expected to be either Space or World.
        _save_bd5(space.as_base(), filename, group_index, object_name, spatial_unit, time_unit, trunc, with_radius)
