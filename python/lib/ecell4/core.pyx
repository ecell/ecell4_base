from libcpp.string cimport string

cdef string tostring(ustr):
    if isinstance(ustr, unicode):
        return <string>(ustr.encode('UTF-8'))
    else:
        return <string>(ustr)

include "RandomNumberGenerator.pxi"
include "UnitSpecies.pxi"
include "Species.pxi"
include "ReactionRule.pxi"
include "Space.pxi"
include "CompartmentSpace.pxi"  #XXX
include "ParticleSpace.pxi"  #XXX
include "Model.pxi"
include "NetworkModel.pxi"
include "NetfreeModel.pxi"
include "Real3.pxi"
include "Integer3.pxi"
include "Particle.pxi"
include "Voxel.pxi"
include "observers.pxi"
include "shapes.pxi"

from cython.operator cimport dereference as deref

cdef shared_ptr[Cpp_Model]* Cpp_Model_from_Model(m):
    if isinstance(m, Model):
        return (<Model>m).thisptr
    elif isinstance(m, NetworkModel):
        return address(<shared_ptr[Cpp_Model]&>(deref((<NetworkModel>m).thisptr)))
    elif isinstance(m, NetfreeModel):
        return address(<shared_ptr[Cpp_Model]&>(deref((<NetfreeModel>m).thisptr)))
    else:
        raise ValueError, ("a wrong argument was given [%s]." % (type(m))
            + " the first argument must be Model, NetworkModel or NetfreeModel")

cimport extras

def load_version_information(string filename):
    return extras.load_version_information(filename)

cimport functions

def cbrt(Real x):
    return functions.cbrt(x)

cimport types

N_A = types.N_A
epsilon = types.epsilon
