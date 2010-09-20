
import math
import numpy
import scipy
import myrandom

import _gfrd

Pi = scipy.pi
Pi2 = scipy.pi * 2.0
PiSqrt = math.sqrt(scipy.pi)

N_A = 6.0221367e23
INF = numpy.inf

ZEROPOS = numpy.array([0., 0., 0.])
NOWHERE = numpy.array((INF, INF, INF))

SAFETY = 1.0 + 1e-5

# Tolerance used for float comparison functions. Oversimplifying: two floats a 
# and b are considered to be equal if abs(a - b) < TOLERANCE * abs(a).
TOLERANCE = 1e-7

# Multiplication factor used for seperating 2 particles or a particle and a 
# surface after unbinding.
MINIMAL_SEPARATION_FACTOR = 1.0 + TOLERANCE
  
# Float comparison functions.
def feq(a, b, typical=1, tolerance=TOLERANCE):
    # Return True if a and b are equal, subject to given tolerances.  
    # Float comparison.

    # Also see numpy.allclose().

    # The (relative) tolerance must be positive and << 1.0

    # Instead of specifying an absolute tolerance, you can speciy a 
    # typical value for a or b. The absolute tolerance is then the 
    # relative tolerance multipied by this typical value, and will be 
    # used when comparing a value to zero. By default, the typical 
    # value is 1.

    return abs(a - b) < tolerance * (typical + min(abs(a), abs(b)))


def fgreater(a, b, typical=1, tolerance=TOLERANCE):
    # Return True if a is greater than b, subject to given tolerances.  
    # Float comparison.

    return a - b > tolerance * (typical + min(abs(a), abs(b)))


def fless(a, b, typical=1, tolerance=TOLERANCE):
    # Return True if a is less than b, subject to given tolerances.  
    # Float comparison.

    return b - a > tolerance * (typical + min(abs(a), abs(b)))


def fgeq(a, b, typical=1, tolerance=TOLERANCE):
    # Return True if a is greater or equal than b, subject to given 
    # tolerances. Float comparison.

    diff = a - b
    barrier = tolerance * (typical + min(abs(a), abs(b)))
    # Try both 'greater than' and equality.
    return diff > barrier or abs(diff) < barrier


def fleq(a, b, typical=1, tolerance=TOLERANCE):
    # Return True if a is less than or equal than b, subject to given 
    # tolerances. Float comparison.

    diff = b - a
    barrier = tolerance * (typical + min(abs(a), abs(b)))
    # Try both 'less than' and equality.
    return diff > barrier or abs(diff) < barrier

def per_M_to_m3(rate):
    """Convert a reaction rate from units 'per molar per second' to 
    units 'meters^3 per second'.

    """
    return rate / (1000 * N_A)

def per_microM_to_m3(rate):
    """Convert a reaction rate from units 'per micromolar per second' to 
    units 'meters^3 per second'.

    """
    return per_M_to_m3(rate * 1e6)

def M_to_per_m3(molar):
    """Convert a concentration from units 'molar' to units 'per 
    meters^3'.

    """
    return molar * (1000 * N_A)

def microM_to_per_m3(micromolar):
    """Convert a concentration from units 'micromolar' to units 'per 
    meters^3'.

    """
    return M_to_per_m3(micromolar / 1e6) 

def mean_arrival_time(r, D):
    return (r * r) / (6.0 * D)

def uniq(l):
    nset = {}
    map(nset.__setitem__, l, [])
    return nset.keys()

cyclic_transpose = _gfrd.cyclic_transpose

def distance_sq_array_simple(position1, positions, fsize = None):
    return numpy.square(positions - position1).sum(1)

def distance_array_simple(position1, positions, fsize = None):
    return numpy.sqrt(distance_sq_array_simple(position1, positions))

distance = _gfrd.distance

distance_cyclic = _gfrd.distance_cyclic

def distance_sq_array_cyclic(position1, positions, fsize):
    diff = numpy.abs(positions - position1)
    diff -= numpy.greater(diff, fsize * 0.5) * fsize # transpose
    return numpy.square(diff).sum(1)

def distance_array_cyclic(position1, positions, fsize = 0):
    return numpy.sqrt(distance_sq_array_cyclic(position1, positions, fsize))

def cartesian_to_spherical(c):
    # x, y, z = c
    r = length(c)
    theta = math.acos(c[2] / r)
    phi = math.atan2(c[1], c[0])
    if phi < 0.0:  # atan2 returns [- PI, PI]
        phi += 2.0 * Pi
    return numpy.array([r, theta, phi])

def spherical_to_cartesian(s):
    #FIXME: it's possible that the below is a source of some bias.
    r, theta, phi = s
    sintheta = math.sin(theta)
    return numpy.array([r * math.cos(phi) * sintheta,
                        r * math.sin(phi) * sintheta,
                        r * math.cos(theta)])

def random_unit_vector_s():
    s = numpy.array([1.0, myrandom.uniform(0, Pi), myrandom.uniform(0, Pi2)])
    return s

def random_unit_vector():
    v = [myrandom.uniform(-1,1), myrandom.uniform(-1,1), myrandom.uniform(-1,1)]
    return _gfrd.normalize(v, 1)

def random_vector(r):
    v = [myrandom.uniform(-1,1), myrandom.uniform(-1,1), myrandom.uniform(-1,1)]
    return _gfrd.normalize(v, r)

def random_vector2D(r):
    # Return a random 2D cartesian vector of length r.

    v = [myrandom.uniform(-1,1), myrandom.uniform(-1,1)]
    # Todo. return _gfrd.normalize(v, r)
    v = numpy.array(v)
    norm = numpy.linalg.norm(v)
    return v * (r / norm)

def length(a):
    return _gfrd.length(a)

def normalize(a, l=1):
    return _gfrd.normalize(a, l)

def vector_angle(a, b):
    cosangle = numpy.dot(a, b) / (length(a) * length(b))
    return math.acos(cosangle)

def vector_angle_against_z_axis(b):
    cosangle = b[2] / length(b)
    return math.acos(cosangle)

def crossproduct(a, b):
    M = numpy.array([[   0.0, - a[2],   a[1]],
                     [  a[2],    0.0, - a[0]],
                     [- a[1],   a[0],    0.0]])
    return numpy.dot(M, b)

def crossproduct_against_z_axis(a):
    return numpy.array([- a[1], a[0], 0.0])

def rotate_vector(v, r, alpha):
    # v: vector to rotate
    # r: normalized rotation axis
    # alpha: rotation angle in radian

    cosalpha = math.cos(alpha)
    sinalpha = math.sin(alpha)
    cosalphac = 1.0 - cosalpha

    M = numpy.array([[cosalpha + cosalphac * r[0] * r[0],
                      cosalphac * r[0] * r[1] - r[2] * sinalpha,
                      cosalphac * r[0] * r[2] + r[1] * sinalpha],
                     [cosalphac * r[0] * r[1] + r[2] * sinalpha,
                      cosalpha + cosalphac * r[1] * r[1],
                      cosalphac * r[1] * r[2] - r[0] * sinalpha],
                     [cosalphac * r[0] * r[2] - r[1] * sinalpha,
                      cosalphac * r[1] * r[2] + r[0] * sinalpha,
                      cosalpha + cosalphac * r[2] * r[2]]])

    return numpy.dot(M,v)

def calculate_pair_CoM(pos1, pos2, D1, D2, world_size):
    return _gfrd.calculate_pair_CoM(pos1, pos2, D1, D2, world_size);

apply_boundary = _gfrd.apply_boundary

def permutate(seq):
    # permutate a sequence and return a list of the permutations
    if not seq:
        return [seq] # is an empty sequence
    else:
        temp = []

        for k in range(len(seq)):
            part = seq[:k] + seq[k+1:]
            for m in permutate(part):
                temp.append(seq[k:k+1] + m)
        return temp

def k_D(Dtot, sigma):
    """Calculate the 'pseudo-'reaction rate (kD) caused by diffusion.
    
    kD is equal to 1 divided by the time it takes for two particles to 
    meet each other by diffusion. It is needed when converting from 
    an intrinsic reaction rate to an overall reaction rates or vice 
    versa.

    Example:
        - A + B -> C.

    Arguments:
        - Dtot:
            the diffusion constant of particle A plus the diffusion 
            constant of particle B. Units: meters^2/second.
        - sigma
            the radius of particle A plus the radius of particle B. 
            Units: meters.

    This function is only available for reaction rules in 3D.

    """
    return 4.0 * numpy.pi * Dtot * sigma

def k_a(kon, kD):
    """Convert an overall reaction rate (kon) for a binding/annihilation 
    reaction rule to an intrinsic reaction rate (ka).

    Example:
        - A + B -> C
            binding reaction rule
        - A + B -> 0
            annihilation reaction rule

    Arguments:
        - kon
            the overall reaction rate for the reaction rule. Units: 
            meters^3/second.
        - kD
            the 'pseudo-'reaction rate caused by the diffusion of 
            particles A and B. See the function k_D(). Units: 
            meters^3/second.

    """
    if kon > kD:
        raise RuntimeError, 'kon > kD.'
    ka = 1. / ((1. / kon) - (1. / kD))
    return ka

def k_d(koff, kon, kD):
    """Convert an overall reaction rate (koff) for an unbinding reaction 
    rule to an intrinsic reaction rate (kd).

    This one is a bit tricky. We consider reaction rules with only 1 
    reactant. In case there is only 1 product also, no conversion in 
    necessary. But when the reaction rule has 2 products, we need to 
    take the reverse reaction rule into account and do the proper 
    conversion.

    Example:
        - C -> A + B
            unbinding reaction rule
        - A + B -> C
            reverse reaction rule

    Arguments:
        - koff
            the overall reaction rate for the unbinding reaction rule.  
            Units: meters^3/second.
        - kon
            the overall reaction rate for the reverse reaction rule. 
            Units: meters^3/second.
        - kD
            the 'pseudo-'reaction rate caused by the diffusion of 
            particles A and B. See the function k_D(). Units: 
            meters^3/second.

    """
    ka = k_a(kon, kD)
    kd = k_d_using_ka(koff, ka, kD)
    return kd

def k_d_using_ka(koff, ka, kD):
    """Convert an overall reaction rate (koff) for an unbinding reaction 
    rule to an intrinsic reaction rate (kd).

    Similar to the function k_d(), but expects an intrinsic rate (ka) 
    instead of an overall rate (kon) for the reversed reaction rule as 
    the second argument.

    """
    kd =  koff * (1 + float(ka) / kD)
    return kd

def k_on(ka, kD):
    """Convert an intrinsic reaction rate (ka) for a binding/annihilation 
    reaction rule to an overall reaction rate (kon).

    The inverse of the function k_a().
    
    Rarely needed.

    """
    kon = 1. / ((1. / kD) + (1. / ka))  # m^3/s
    return kon

def k_off(kd, kon, kD):
    """Convert an intrinsic reaction rate (kd) for an unbinding reaction 
    rule to an overall reaction rate (koff).

    The inverse of the function k_d().

    Rarely needed.

    """
    ka = k_a(kon, kD) 
    koff = k_off_using_ka(kd, ka, kD)
    return koff

def k_off_using_ka(kd, ka, kD):
    """Convert an intrinsic reaction rate (kd) for an unbinding reaction 
    rule to an overall reaction rate (koff).

    Similar to the function k_off(), but expects an intrinsic rate 
    (ka) instead of an overall rate (kon) as the second argument.

    Rarely needed.

    """
    koff = 1. / (float(ka) / (kd * kD) + (1. / kd))
    return koff

def C2N(c, V):
    return c * V * N_A  # round() here?

