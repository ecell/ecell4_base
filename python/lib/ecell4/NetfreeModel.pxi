from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address

from libcpp.vector cimport vector
from libcpp.map cimport map
from shared_ptr cimport shared_ptr


cdef class NetfreeModel:
    """A netfree model class.

    NetfreeModel()

    """

    def __init__(self):
        """Constructor."""
        pass

    def __cinit__(self):
        # self.thisptr = new NetfreeModel()
        self.thisptr = new shared_ptr[Cpp_NetfreeModel](new Cpp_NetfreeModel())

    def __dealloc__(self):
        del self.thisptr

    def add_species_attribute(self, Species sp):
        """add_species_attribute(sp)

        Add a species attribute to the bottom.

        Parameters
        ----------
        sp : Species
            A new species with attributes.

        """
        self.thisptr.get().add_species_attribute(deref(sp.thisptr))

    def has_species_attribute(self, Species sp):
        """has_species_attribute(sp) -> bool

        Return if the given species can be attributed or not.

        """
        return self.thisptr.get().has_species_attribute(deref(sp.thisptr))

    def remove_species_attribute(self, Species sp):
        """remove_species_attribute(sp)

        Remove the species attribute.

        """
        self.thisptr.get().remove_species_attribute(deref(sp.thisptr))

    def add_reaction_rule(self, ReactionRule rr):
        """add_reaction_rule(rr)

        Add a new reaction rule.

        Parameters
        ----------
        rr : ReactionRule
            A new reaction rule.

        """
        self.thisptr.get().add_reaction_rule(deref(rr.thisptr))

    def remove_reaction_rule(self, ReactionRule rr):
        """remove_reaction_rule(rr)

        Remove a reaction rule.

        """
        self.thisptr.get().remove_reaction_rule(deref(rr.thisptr))

    def has_reaction_rule(self, ReactionRule rr):
        """has_reaction_rule(rr) -> bool

        Return if the given reaction rule is existing or not.

        """
        return self.thisptr.get().has_reaction_rule(deref(rr.thisptr))

    def num_reaction_rules(self):
        """Return a number of reaction rules contained in the model."""
        return self.thisptr.get().num_reaction_rules()

    def apply_species_attributes(self, Species sp):
        """apply_species_attributes(sp) -> Species

        Return a species with attributes.

        Parameters
        ----------
        sp : Species
            An original species.

        Returns
        -------
        Species:
            A new species attributed by species attributes in the model.

        """
        cdef Cpp_Species retval = self.thisptr.get().apply_species_attributes(
            deref(sp.thisptr))
        return Species_from_Cpp_Species(address(retval))

    # def create_species(self, string name):
    #     cdef Cpp_Species retval = self.thisptr.get().create_species(name)
    #     return Species_from_Cpp_Species(address(retval))

    def reaction_rules(self):
        """Return a list of reaction rules contained in the model."""
        cdef vector[Cpp_ReactionRule] c_rr_vector = self.thisptr.get().reaction_rules()
        retval = []
        cdef vector[Cpp_ReactionRule].iterator it = c_rr_vector.begin()
        while it != c_rr_vector.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(
                <Cpp_ReactionRule*>(address(deref(it)))))
            inc(it)
        return retval

    def species_attributes(self):
        """Return a list of species attributes contained in the model."""
        cdef vector[Cpp_Species] species = self.thisptr.get().species_attributes()
        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(Species_from_Cpp_Species(
                <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def list_species(self):
        """Return a list of species, contained in reaction rules in the model."""
        cdef vector[Cpp_Species] species = self.thisptr.get().list_species()
        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(Species_from_Cpp_Species(
                <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def query_reaction_rules(self, Species sp1, Species sp2 = None):
        """query_reaction_rules(sp1, sp2=None) -> [ReactionRule]

        Query and return a list of reaction rules, which have the given species
        as their reactants.

        Parameters
        ----------
        sp1 : Species
            The first reactant
        sp2 : Species
            The second reactant. This is for querying second order reaction rules.

        Returns
        -------
        list:
            A list of ``ReactionRule``s.

        """
        cdef vector[Cpp_ReactionRule] rules
        if sp2 is None:
            rules = self.thisptr.get().query_reaction_rules(
                deref(sp1.thisptr))
        else:
            rules = self.thisptr.get().query_reaction_rules(
                deref(sp1.thisptr), deref(sp2.thisptr))
        retval = []
        cdef vector[Cpp_ReactionRule].iterator it = rules.begin()
        while it != rules.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(
                <Cpp_ReactionRule*>(address(deref(it)))))
            inc(it)
        return retval

    def add_species_attributes(self, attrs):
        """add_species_attributes(attrs)

        Extend a list of species attributes to the bottom.

        Parameters
        ----------
        attrs : list
            A list of new ``Species`` with attributes.

        """
        cdef vector[Cpp_Species] species
        for sp in attrs:
            species.push_back(deref((<Species>sp).thisptr))
        self.thisptr.get().add_species_attributes(species)

    def add_reaction_rules(self, rrs):
        """add_reaction_rules(rrs)

        Add a list of new reaction rules.

        Parameters
        ----------
        rrs : list
            A list of new ``ReactionRule``s.

        """
        cdef vector[Cpp_ReactionRule] reaction_rules
        for rr in rrs:
            reaction_rules.push_back(deref((<ReactionRule>rr).thisptr))
        self.thisptr.get().add_reaction_rules(reaction_rules)

    def expand(self, seeds, max_itr=None, max_stoich=None):
        """expand(seeds, max_itr=None, max_stoich=None) -> Model

        Expand a rule-based model into a network model.

        Parameters
        ----------
        seeds : list
            A list of ``Species`` which gives seeds.
        max_itr : Integer
            A maximum number of iterations to generate new products.
        max_stoich : Integer
            A maximum stoichiometry of ``UnitSpecies`` in a ``Species``.

        Returns
        -------
        Model:
            A network model.

        """
        cdef vector[Cpp_Species] _seeds
        cdef map[Cpp_Species, Integer] _max_stoich
        for sp in seeds:
            if not isinstance(sp, Species):
                raise ValueError(
                    'seeds must be given as a list of Species.'
                    + ' {0} given.'.format(repr(sp)))
            _seeds.push_back(deref((<Species>sp).thisptr))

        if max_stoich is not None:
            for sp, n in max_stoich.items():
                _max_stoich[deref((<Species>sp).thisptr)] = <Integer>n
            return Model_from_Cpp_Model(
                self.thisptr.get().expand(_seeds, <Integer>max_itr, _max_stoich))
        elif max_itr is not None:
            return Model_from_Cpp_Model(
                self.thisptr.get().expand(_seeds, <Integer>max_itr))
        else:
            return Model_from_Cpp_Model(
                self.thisptr.get().expand(_seeds))

    def add_parameter(self, Species sp):
        """add_parameter(sp)

        This is for the tentative implementation of parameters.
        This might be deprecated.

        """
        self.thisptr.get().add_parameter(deref(sp.thisptr))

    def add_parameters(self, attrs):
        """add_parameters(attrs)

        This is for the tentative implementation of parameters.
        This might be deprecated.

        """
        cdef vector[Cpp_Species] species
        for sp in attrs:
            species.push_back(deref((<Species>sp).thisptr))
        self.thisptr.get().add_parameters(species)

    def parameters(self):
        """parameters()

        This is for the tentative implementation of parameters.
        This might be deprecated.

        """
        cdef vector[Cpp_Species] species = self.thisptr.get().parameters()
        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(Species_from_Cpp_Species(
                <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

cdef NetfreeModel NetfreeModel_from_Cpp_NetfreeModel(
    shared_ptr[Cpp_NetfreeModel] m):
    r = NetfreeModel()
    r.thisptr.swap(m)
    return r
