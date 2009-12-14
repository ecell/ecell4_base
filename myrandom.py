from _gfrd import RandomNumberGenerator
import numpy

__all__ = (
    'shuffle',
    'uniform',
    'normal',
    'seed',
    )

rng = RandomNumberGenerator()

def uniform(min=0.0, max=1.0, size=None):
    global rng
    if size is None:
        return rng.uniform(min, max)
    retval = numpy.ndarray(dtype=numpy.float64, shape=size)
    f = retval.flat
    for i in range(len(f)):
        f[i] = rng.uniform(min, max)
    return retval

normal = rng.normal
seed = rng.seed

def shuffle(seq):
    for i in reversed(range(0, len(seq))):
        j = rng.uniform_int(0, i)
        seq[i], seq[j] = seq[j], seq[i]

