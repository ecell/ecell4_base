from cython.operator cimport dereference as deref
from libcpp.string cimport string
cimport util


cdef class UnitSpecies:
    """A class representing an unit of species.

    UnitSpecies(name=None)

    See Also
    --------
    Species

    """

    def __init__(self, name=None):
        """Constructor.

        Parameters
        ----------
        name : str, optional
            A name.

        """
        pass

    def __cinit__(self, name=None):
        if name is None:
            self.thisptr = new Cpp_UnitSpecies()
        else:
            self.thisptr = new Cpp_UnitSpecies(tostring(name))

    def __dealloc__(self):
        del self.thisptr

    def __richcmp__(UnitSpecies self, UnitSpecies rhs, int op):
        cdef int compare
        if deref(self.thisptr) > deref(rhs.thisptr):
            compare = 1
        elif deref(self.thisptr) < deref(rhs.thisptr):
            compare = -1
        else: # self == rhs
            compare = 0
        return util.richcmp_helper(compare, op)

    def __hash__(self):
        return hash(self.thisptr.serial())

    def serial(self):
        """Return the serial, which consists of a name and sites."""
        return self.thisptr.serial().decode('UTF-8')

    def name(self):
        """Return a name."""
        return self.thisptr.name().decode('UTF-8')

    def add_site(self, name, state, bond):
        """add_site(name, state, bond)

        Add a new site.

        Parameters
        ----------
        name : str
            A name of the site
        state : str
            A state name of the site
        bond : str
            A bond of the site.

        """
        return self.thisptr.add_site(tostring(name), tostring(state), tostring(bond))

    def deserialize(self, serial):
        """deserialize(serial)

        Deserialize the given serial, and load it.

        Parameters
        ----------
        serial : str
            A serial

        """
        self.thisptr.deserialize(tostring(serial))

cdef UnitSpecies UnitSpecies_from_Cpp_UnitSpecies(Cpp_UnitSpecies *sp):
    cdef Cpp_UnitSpecies *new_obj = new Cpp_UnitSpecies(deref(sp))
    r = UnitSpecies('')
    del r.thisptr
    r.thisptr = new_obj
    return r
