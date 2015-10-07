from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.string cimport string
from cython cimport address
cimport util

cimport context


cdef class Species:
    """A class representing a type of molecules with attributes.

    Species(serial=None, radius=None, D=None, location=None)

    """

    def __init__(self, serial=None, radius=None, D=None, location=None):
        """Constructor.

        Args:
            serial (str, optional): The serial name.
            radius (str, optional): The radius of a molecule.
                This must be given as a string.
            D (str, optional): The diffusion rate of a molecule.
                This must be given as a string.
            location (str, optional): The location of a molecule.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, serial=None, radius=None, D=None, location=None):
        if serial is None:
            self.thisptr = new Cpp_Species()
        elif radius is not None and D is not None:
            if location is None:
                self.thisptr = new Cpp_Species(
                    tostring(serial),
                    tostring(radius),
                    tostring(D))
            else:
                self.thisptr = new Cpp_Species(
                    tostring(serial),
                    tostring(radius),
                    tostring(D),
                    tostring(location))
        else:
            self.thisptr = new Cpp_Species(tostring(serial)) #XXX:

    def __dealloc__(self):
        del self.thisptr

    def __richcmp__(Species self, Species rhs, int op):
        cdef int compare
        if deref(self.thisptr) > deref(rhs.thisptr):
            compare = 1
        elif deref(self.thisptr) < deref(rhs.thisptr):
            compare = -1
        else: # self == rhs
            compare = 0
        return util.richcmp_helper(compare, op)

    def __hash__(self):
        return hash(self.thisptr.serial().decode('UTF-8'))

    def serial(self):
        """Return the serial name as an unicode string."""
        return self.thisptr.serial().decode('UTF-8')

    def get_attribute(self, name):
        """get_attribute(name) -> str

        Return an attribute as an unicode string.
        If no corresponding attribute is found, raise an error.

        Args:
            name (str): The name of an attribute.

        Returns:
            str: The value of the attribute.

        """
        return self.thisptr.get_attribute(
            tostring(name)).decode('UTF-8')

    def set_attribute(self, name, value):
        """set_attribute(name, value)

        Set an attribute.
        If existing already, the attribute will be overwritten.

        Args:
            name (str): The name of an attribute.
            value (str): The value of an attribute.

        """
        self.thisptr.set_attribute(tostring(name), tostring(value))

    def remove_attribute(self, name):
        """remove_attribute(name)

        Remove an attribute.
        If no corresponding attribute is found, raise an error.

        Args:
            name (str): The name of an attribute to be removed.

        """
        self.thisptr.remove_attribute(tostring(name))

    def has_attribute(self, name):
        """has_attribute(name) -> bool

        Return if the attribute exists or not.

        Args:
            name (str): The name of an attribute.

        Returns:
            bool: True if the attribute exists, False otherwise.

        """
        return self.thisptr.has_attribute(tostring(name))

    def list_attributes(self):
        """list_attributes() -> [(str, str)]

        List all attributes.

        Returns:
            list: A list of pairs of name and value.
                ``name`` and ``value`` are given as unicode strings.

        """
        retval = self.thisptr.list_attributes()
        return [(key.decode('UTF-8'), value.decode('UTF-8'))
            for key, value in retval]

    def add_unit(self, UnitSpecies usp):
        """add_unit(usp)

        Append an ``UnitSpecies`` to the end.

        Args:
            usp (UnitSpecies): An ``UnitSpecies`` to be added.

        """
        self.thisptr.add_unit(deref(usp.thisptr))

    def count(self, Species pttrn):
        """count(pttrn) -> Integer

        Count the number of matches for a pattern given as a ``Species``.

        Args:
            pttrn (Species): A pattern to be count.

        Returns:
            Integer: The number of matches.

        """
        return self.thisptr.count(deref(pttrn.thisptr))

    def units(self):
        """units() -> [UnitSpecies]

        Return a list of all ``UnitSpecies`` contained.

        """
        cdef vector[Cpp_UnitSpecies] usps = self.thisptr.units()
        retval = []
        cdef vector[Cpp_UnitSpecies].iterator it = usps.begin()
        while it != usps.end():
            retval.append(UnitSpecies_from_Cpp_UnitSpecies(
            <Cpp_UnitSpecies*>(address(deref(it)))))
            inc(it)
        return retval

    def num_units(self):
        """num_units() -> Integer

        Return the number of ``UnitSpecies``.

        """
        return self.thisptr.num_units()

    def deserialize(self, serial):
        """deserialize(serial)

        Reset the serial. All attributes will be kept.

        Args:
            serial (str): A new serial as an unicode string.

        """
        self.thisptr.deserialize(tostring(serial))

cdef Species Species_from_Cpp_Species(Cpp_Species *sp):
    cdef Cpp_Species *new_obj = new Cpp_Species(deref(sp))
    r = Species()
    del r.thisptr
    r.thisptr = new_obj
    return r

def spmatch(Species pttrn, Species sp):
    """spmatch(pttrn, sp) -> bool

    Return if a pattern matches the target ``Species`` or not.

    Args:
        pttrn (Species): A pattern.
        sp (Species): A target.

    Return:
        bool: True if ``pttrn`` matches ``sp`` at least one time, False otherwise.

    """
    return context.spmatch(deref(pttrn.thisptr), deref(sp.thisptr))

def count_spmatches(Species pttrn, Species sp):
    """count_spmatches(pttrn, sp) -> Integer

    Count the number of matches for a pattern given as a ``Species``.

    Args:
        pttrn (Species): A pattern.
        sp (Species): A target.

    Return:
        Integer: The number of matches.

    Note:
        Rather use ``Species.count``.

    """
    return context.count_spmatches(deref(pttrn.thisptr), deref(sp.thisptr))

def format_species(Species sp):
    """format_species(sp) -> Species

    Return a species uniquely reformatted.

    """
    cdef Cpp_Species newsp = context.format_species(deref(sp.thisptr))
    return Species_from_Cpp_Species(address(newsp))

def unique_serial(Species sp):
    """unique_serial(sp) -> str

    Return a serial of a species uniquely reformatted.
    This is equivalent to call ``format_species(sp).serial()``

    """
    return context.unique_serial(deref(sp.thisptr)).decode('UTF-8')
