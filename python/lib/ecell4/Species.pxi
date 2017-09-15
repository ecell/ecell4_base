from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.string cimport string
from cython cimport address

cimport util
cimport context

import numbers
from cpython cimport bool as bool_t

cdef boost_get_from_Cpp_Species_value_type(Cpp_Species_value_type value):
    cdef string* value_str = boost_get[string, string, Real, Integer, bool](address(value))
    if value_str != NULL:
        return deref(value_str).decode('UTF-8')
    cdef Real* value_real = boost_get[Real, string, Real, Integer, bool](address(value))
    if value_real != NULL:
        return deref(value_real)
    cdef Integer* value_int = boost_get[Integer, string, Real, Integer, bool](address(value))
    if value_int != NULL:
        return deref(value_int)
    cdef bool* value_bool = boost_get[bool, string, Real, Integer, bool](address(value))
    if value_bool != NULL:
        return deref(value_bool)

    raise RuntimeError('Never get here. Unsupported return type was given.')

cdef class Species:
    """A class representing a type of molecules with attributes.

    Species(serial=None, radius=None, D=None, location=None)

    """

    def __init__(self, serial=None, radius=None, D=None, location=None):
        """Constructor.

        Parameters
        ----------
        serial : str, optional
            The serial name.
        radius : float, optional
            The radius of a molecule.
        D : foat, optional
            The diffusion rate of a molecule.
        location : str, optional
            The location of a molecule.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, serial=None, radius=None, D=None, location=None):
        if serial is None:
            self.thisptr = new Cpp_Species()
        elif radius is None:
            self.thisptr = new Cpp_Species(tostring(serial)) #XXX:
        elif D is None:
            raise ValueError(
                'D must be given. D is not optional when radius is given.')
        elif location is None:
            if isinstance(radius, str) and isinstance(D, str):
                self.thisptr = new Cpp_Species(
                    tostring(serial), tostring(radius), tostring(D))
            elif isinstance(radius, numbers.Real) and isinstance(D, numbers.Real):
                self.thisptr = new Cpp_Species(
                    tostring(serial), <Real>radius, <Real>D)
            else:
                raise TypeError('radius and D must be float.')
        else:
            if isinstance(radius, str) and isinstance(D, str):
                self.thisptr = new Cpp_Species(
                    tostring(serial), tostring(radius), tostring(D), tostring(location))
            elif isinstance(radius, numbers.Real) and isinstance(D, numbers.Real):
                self.thisptr = new Cpp_Species(
                    tostring(serial), <Real>radius, <Real>D, tostring(location))
            else:
                raise TypeError('radius and D must be float.')

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

    # def get_attribute(self, name):
    #     """get_attribute(name) -> str

    #     Return an attribute as an unicode string.
    #     If no corresponding attribute is found, raise an error.

    #     Parameters
    #     ----------
    #     name : str
    #         The name of an attribute.

    #     Returns
    #     -------
    #     value : str
    #         The value of the attribute.

    #     """
    #     return self.thisptr.get_attribute(
    #         tostring(name)).decode('UTF-8')

    def get_attribute(self, name):
        """get_attribute(name) -> str, float, int, or bool

        Return an attribute.
        If no corresponding attribute is found, raise an error.

        Parameters
        ----------
        name : str
            The name of an attribute.

        Returns
        -------
        value : str, float, int, or bool
            The value of the attribute.

        """
        return boost_get_from_Cpp_Species_value_type(self.thisptr.get_attribute(tostring(name)))

    def set_attribute(self, name, value):
        """set_attribute(name, value)

        Set an attribute.
        If existing already, the attribute will be overwritten.

        Parameters
        ----------
        name : str
            The name of an attribute.
        value : str, float, int, or bool
            The value of an attribute.

        """
        if isinstance(value, str):
            self.thisptr.set_attribute(tostring(name), tostring(value))
        elif isinstance(value, bool_t):
            self.thisptr.set_attribute(tostring(name), <bool> value)
        elif isinstance(value, numbers.Integral):
            self.thisptr.set_attribute(tostring(name), <Integer> value)
        elif isinstance(value, numbers.Real):
            self.thisptr.set_attribute(tostring(name), <Real> value)
        else:
            raise TypeError(
                'Type [{}] is not supported. str, int, float or bool must be given.'.format(
                    type(value)))

    def remove_attribute(self, name):
        """remove_attribute(name)

        Remove an attribute.
        If no corresponding attribute is found, raise an error.

        Parameters
        ----------
        name : str
            The name of an attribute to be removed.

        """
        self.thisptr.remove_attribute(tostring(name))

    def has_attribute(self, name):
        """has_attribute(name) -> bool

        Return if the attribute exists or not.

        Parameters
        ----------
        name : str
            The name of an attribute.

        Returns
        -------
        bool:
            True if the attribute exists, False otherwise.

        """
        return self.thisptr.has_attribute(tostring(name))

    def list_attributes(self):
        """list_attributes() -> [(str, str)]

        List all attributes.

        Returns
        -------
        list:
            A list of pairs of name and value.
            ``name`` and ``value`` are given as unicode strings.

        """
        cdef vector[pair[string, Cpp_Species_value_type]] attrs = self.thisptr.list_attributes()
        res = []
        cdef vector[pair[string, Cpp_Species_value_type]].iterator it = attrs.begin()
        while it != attrs.end():
            res.append((deref(it).first.decode('UTF-8'), boost_get_from_Cpp_Species_value_type(deref(it).second)))
            inc(it)
        return res

    def add_unit(self, UnitSpecies usp):
        """add_unit(usp)

        Append an ``UnitSpecies`` to the end.

        Parameters
        ----------
        usp : UnitSpecies
            An ``UnitSpecies`` to be added.

        """
        self.thisptr.add_unit(deref(usp.thisptr))

    def count(self, Species sp):
        """count(sp) -> Integer

        Count the number of matches for a target given as a ``Species``.

        Parameters
        ----------
        sp : Species
            A target to be count.

        Returns
        -------
        Integer:
            The number of matches.

        """
        return self.thisptr.count(deref(sp.thisptr))

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

    def D(self, value):
        """D(string) -> Species

        set attribute 'D', and return self.

        """
        cdef Cpp_Species *sp = self.thisptr.D_ptr(tostring(value))
        assert sp == self.thisptr
        return self

    def radius(self, value):
        """radius(string) -> Species

        set attribute 'radius', and return self.

        """
        cdef Cpp_Species *sp = self.thisptr.radius_ptr(tostring(value))
        assert sp == self.thisptr
        return self

    def location(self, value):
        """location(string) -> Species

        set attribute 'location', and return self.

        """
        cdef Cpp_Species *sp = self.thisptr.location_ptr(tostring(value))
        assert sp == self.thisptr
        return self

    def __reduce__(self):
        return (__rebuild_species, (self.serial(), self.list_attributes()))

def __rebuild_species(serial, attrs):
    sp = Species(serial)
    for key, val in attrs:
        sp.set_attribute(key, val)
    return sp

cdef Species Species_from_Cpp_Species(Cpp_Species *sp):
    cdef Cpp_Species *new_obj = new Cpp_Species(deref(sp))
    r = Species()
    del r.thisptr
    r.thisptr = new_obj
    return r

def spmatch(Species pttrn, Species sp):
    """spmatch(pttrn, sp) -> bool

    Return if a pattern matches the target ``Species`` or not.

    Parameters
    ----------
    pttrn : Species
        A pattern.
    sp : Species
        A target.

    Returns
    -------
    bool:
        True if ``pttrn`` matches ``sp`` at least one time, False otherwise.

    """
    return context.spmatch(deref(pttrn.thisptr), deref(sp.thisptr))

def count_spmatches(Species pttrn, Species sp):
    """count_spmatches(pttrn, sp) -> Integer

    Count the number of matches for a pattern given as a ``Species``.

    Parameters
    ----------
    pttrn : Species
        A pattern.
    sp : Species
        A target.

    Returns
    -------
    Integer:
        The number of matches.

    Notes
    -----
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
