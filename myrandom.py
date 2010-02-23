from _gfrd import RandomNumberGenerator
import numpy

__all__ = (
    'shuffle',
    'uniform',
    'normal',
    'seed',
    )

rng = RandomNumberGenerator.create()

def uniform(min=0.0, max=1.0, size=None):
     global rng
     return rng.uniform(min, max)

normal = rng.normal
seed = rng.seed

def shuffle(seq):
    for i in reversed(range(0, len(seq))):
        j = rng.uniform_int(0, i)
        seq[i], seq[j] = seq[j], seq[i]

def choice(a, b):
    '''Return a or b with 50% probability each.

    '''
    return uniform() > 0.5 and a or b
