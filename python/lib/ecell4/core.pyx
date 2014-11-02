include "UnitSpecies.pxi"
include "Species.pxi"
include "Ratelaw.pxi"
include "ReactionRule.pxi"
include "Space.pxi"
include "CompartmentSpace.pxi"
include "ParticleSpace.pxi"
include "Model.pxi"
include "NetworkModel.pxi"
include "NetfreeModel.pxi"
include "Real3.pxi"
include "Integer3.pxi"
include "Particle.pxi"
include "RandomNumberGenerator.pxi"
include "Voxel.pxi"
include "Observer.pxi"
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
