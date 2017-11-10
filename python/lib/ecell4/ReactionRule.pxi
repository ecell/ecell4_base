from cython.operator cimport dereference as deref, preincrement as inc
from cython cimport address

cimport create_reaction_rule as crr


class ReactionRulePolicy(object):
    """A wrapper of ReactionRule::policy_type"""

    def __init__(self, value):
        self.__value = value

    def get(self):
        return self.__value

cdef class ReactionRule:
    """A class representing a reaction rule between ``Species``.

    ReactionRule(reactants=None, products=None, k=None)

    """

    # STRICT = Cpp_STRICT
    # IMPLICIT = Cpp_IMPLICIT
    # DESTROY = Cpp_DESTROY
    STRICT = ReactionRulePolicy(Cpp_STRICT)
    IMPLICIT = ReactionRulePolicy(Cpp_IMPLICIT)
    DESTROY = ReactionRulePolicy(Cpp_DESTROY)

    def __init__(self, reactants=None, products=None, k=None):
        """Constructor.

        Parameters
        ----------
        reactants : list, optional
            A list of reactant ``Species``.
        products : list, optional
            A list of product ``Species``.
        k : float, optional
            A kinetic rate constant.

        """
        pass  # XXX: Only used for doc string

    def __cinit__(self, reactants=None, products=None, k=None):
        cdef vector[Cpp_Species] cpp_reactants
        cdef vector[Cpp_Species] cpp_products

        if products is None:
            self.thisptr = new Cpp_ReactionRule()
        else:
            for sp in reactants:
                cpp_reactants.push_back(deref((<Species>sp).thisptr))
            for sp in products:
                cpp_products.push_back(deref((<Species>sp).thisptr))

            if k is None:
                self.thisptr = new Cpp_ReactionRule(cpp_reactants, cpp_products)
            else:
                self.thisptr = new Cpp_ReactionRule(cpp_reactants, cpp_products, k)

    def __dealloc__(self):
        del self.thisptr

    def k(self):
        """Return the kinetic rate constant as a float value."""
        return self.thisptr.k()

    def set_k(self, Real k):
        """set_k(k)

        Set a kinetic rate constant.

        Parameters
        ----------
        k : float
            A kinetic rate constant.

        """
        self.thisptr.set_k(k)

    def reactants(self):
        """List all reactants.

        Returns
        -------
        list:
            A list of reactant ``Species``.

        """
        cdef vector[Cpp_Species] reactants = self.thisptr.reactants()
        retval = []
        cdef vector[Cpp_Species].iterator it = reactants.begin()
        while it != reactants.end():
            retval.append(
                Species_from_Cpp_Species(<Cpp_Species*>address(deref(it))))
            inc(it)
        return retval

    def products(self):
        """List all products.

        Returns
        -------
        list:
            A list of product ``Species``.

        """
        cdef vector[Cpp_Species] products = self.thisptr.products()
        retval = []
        cdef vector[Cpp_Species].iterator it = products.begin()
        while it != products.end():
            retval.append(
                Species_from_Cpp_Species(<Cpp_Species*>address(deref(it))))
            inc(it)
        return retval

    def add_reactant(self, Species sp):
        """add_reactant(sp)

        Append a reactant to the end.

        Parameters
        ----------
        sp : Species
            A new reactant.

        """
        self.thisptr.add_reactant(deref(sp.thisptr))

    def add_product(self, Species sp):
        """add_product(sp)

        Append a product to the end.

        Parameters
        ----------
        sp : Species
            A new product.

        """
        self.thisptr.add_product(deref(sp.thisptr))

    def as_string(self):
        """as_string() -> str

        Return an unicode string describing this object.

        Returns
        -------
        str:
            An unicode string describing this object.

        Examples
        --------
        The string consists of a list of reactants, a list of products,
        and a kinetic rate constant.

        >>> rr = ReactionRule([Species("A"), Species("B")], [Species("C")], 1.0)
        >>> rr.as_string()
        u'A+B>C|1'
        """
        return self.thisptr.as_string().decode('UTF-8')

    def policy(self):
        """policy() -> int

        Return a policy for the rule-based modeling.

        """
        return self.thisptr.policy()

    def set_policy(self, policy):
        """set_policy(policy)

        Set a policy for the rule-based modeling.

        Examples
        --------

        >>> rr = ReactionRule()
        >>> rr.set_policy(ReactionRule.STRICT | ReactionRule.DESTROY)

        """
        self.thisptr.set_policy(policy)

    def count(self, reactants):
        """count(reactants) -> Integer

        Count the number of matches for reactants.

        Parameters
        ----------
        reactants : list
            A list of ``Species``. The order of ``reactants``
            is respected.

        Returns
        -------
        Integer:
            The number of matches.

        """
        cdef vector[Cpp_Species] cpp_reactants
        for sp in reactants:
            cpp_reactants.push_back(deref((<Species> sp).thisptr))
        return self.thisptr.count(cpp_reactants)

    def generate(self, reactants):
        """generate(reactants) -> [ReactionRule]

        Generate ``ReactionRule``s from given reactants.

        Parameters
        ----------
        reactants : list
            A list of ``Species``. The order of ``reactants`` is respected.

        Returns
        -------
        list:
            A list of ``ReactionRule``s. The reactants of each
            ``ReactionRule`` are equal to the given ``reactants``.
            If the ``ReactionRule`` does not match the ``reactants``,
            return an empty list.

        Examples
        --------

        >>> rr = ReactionRule([Species("_(b=x)")], [Species("_(b=y)")], 1.0)
        >>> reactants = [Species("A(a^1,b=x).B(a^1,b=x)")]
        >>> [r.as_string() for r in rr.generate(reactants)]
        [u'A(a^1,b=x).B(a^1,b=x)>A(a^1,b=y).B(a^1,b=x)|1',
         u'A(a^1,b=x).B(a^1,b=x)>A(a^1,b=x).B(a^1,b=y)|1']

        """
        cdef vector[Cpp_Species] cpp_reactants
        for sp in reactants:
            cpp_reactants.push_back(deref((<Species> sp).thisptr))
        cdef vector[Cpp_ReactionRule] cpp_rules = self.thisptr.generate(cpp_reactants)
        cdef vector[Cpp_ReactionRule].iterator it1 = cpp_rules.begin()
        retval = []
        while it1 != cpp_rules.end():
            retval.append(ReactionRule_from_Cpp_ReactionRule(address(deref(it1))))
            inc(it1)
        return retval

    def __reduce__(self):
        return (ReactionRule, (self.reactants(), self.products(), self.k()))

cdef ReactionRule ReactionRule_from_Cpp_ReactionRule(Cpp_ReactionRule *rr):
    cdef Cpp_ReactionRule *new_obj = new Cpp_ReactionRule(deref(rr))
    r = ReactionRule()
    del r.thisptr
    r.thisptr = new_obj
    return r

def create_degradation_reaction_rule(Species reactant1, Real k):
    """create_degradation_reaction_rule(reactant1, k) -> ReactionRule

    Create a degradation ``ReactionRule``.

    Parameters
    ----------
    reactant1 : Species
        A reactant to be degradated.
    k : float
        A kinetic parameter.

    Notes
    -----
    This is equivalent to ``ReactionRule([reactant1], [], k)``.

    """
    cdef Cpp_ReactionRule rr = crr.create_degradation_reaction_rule(
        deref(reactant1.thisptr), k)
    return ReactionRule_from_Cpp_ReactionRule(address(rr))

def create_synthesis_reaction_rule(Species product1, Real k):
    """create_synthesis_reaction_rule(product1, k) -> ReactionRule

    Create a synthesis ``ReactionRule``.

    Parameters
    ----------
    product1 : Species
        A product to be synthesized.
    k : float
        A kinetic parameter.

    Notes
    -----
    This is equivalent to ``ReactionRule([], [product1], k)``.

    """
    cdef Cpp_ReactionRule rr = crr.create_synthesis_reaction_rule(
        deref(product1.thisptr), k)
    return ReactionRule_from_Cpp_ReactionRule(address(rr))

def create_unimolecular_reaction_rule(Species reactant1, Species product1, Real k):
    """create_unimolecular_reaction_rule(reactant1, product1, k) -> ReactionRule

    Create an unimolecular ``ReactionRule``.

    Parameters
    ----------
    reactant1 : Species
        A reactant to be modified.
    product1 : Species
        A product.
    k : float
        A kinetic parameter.

    Notes
    -----
    This is equivalent to ``ReactionRule([reactant1], [product1], k)``.

    """
    cdef Cpp_ReactionRule rr = crr.create_unimolecular_reaction_rule(
        deref(reactant1.thisptr), deref(product1.thisptr), k)
    return ReactionRule_from_Cpp_ReactionRule(address(rr))

def create_binding_reaction_rule(
    Species reactant1, Species reactant2, Species product1, Real k):
    """create_binding_reaction_rule(reactant1, reactant2, product1, k) -> ReactionRule

    Create a binding ``ReactionRule``.

    Parameters
    ----------
    reactant1 : Species
        One of two reactants.
    reactant2 : Species
        One of two reactants.
    product1 : Species
        A product.
    k : float
        A kinetic parameter.

    Notes
    -----
    This is equivalent to ``ReactionRule([reactant1, reactant2], [product1], k)``.

    """
    cdef Cpp_ReactionRule rr = crr.create_binding_reaction_rule(
        deref(reactant1.thisptr), deref(reactant2.thisptr),
        deref(product1.thisptr), k)
    return ReactionRule_from_Cpp_ReactionRule(address(rr))

def create_unbinding_reaction_rule(
    Species reactant1, Species product1, Species product2, Real k):
    """create_unbinding_reaction_rule(reactant1, product1, product2, k) -> ReactionRule

    Create an unbinding ``ReactionRule``.

    Parameters
    ----------
    reactant1 : Species
        A reactant.
    product1 : Species
        One of two products.
    product2 : Species
        One of two products.
    k : float
        A kinetic parameter.

    Notes
    -----
    This is equivalent to ``ReactionRule([reactant1], [product1, product2], k)``.

    """
    cdef Cpp_ReactionRule rr = crr.create_unbinding_reaction_rule(
        deref(reactant1.thisptr),
        deref(product1.thisptr), deref(product2.thisptr), k)
    return ReactionRule_from_Cpp_ReactionRule(address(rr))

# def rrmatch(ReactionRule pttrn, reactants):
#     """rrmatch(pttrn, reactants) -> bool
# 
#     Return if a pattern matches the reactants or not.
# 
#     Parameters
#     ----------
#     pttrn : ReactionRule
#         A pattern.
#     reactants : list
#         A list of reactants, ``Species``. The order of reactants is respected.
# 
#     Returns
#     -------
#     bool:
#         True if ``pttrn`` matches ``reactants`` at least one time,
#         False otherwise.
# 
#     """
#     cdef vector[Cpp_Species] cpp_reactants
#     for sp in reactants:
#         cpp_reactants.push_back(deref((<Species> sp).thisptr))
#     return context.rrmatch(deref(pttrn.thisptr), cpp_reactants)
# 
# def count_rrmatches(ReactionRule pttrn, reactants):
#     """count_rrmatches(pttrn, reactants) -> Integer
# 
#     Count the number of matches for a pattern given as a ``ReactionRule``.
# 
#     Parameters
#     ----------
#     pttrn : ReactionRule
#         A pattern.
#     reactants : list
#         A list of reactants, ``Species``. The order of reactants is respected.
# 
#     Returns
#     -------
#     Integer:
#         The number of matches.
# 
#     """
#     cdef vector[Cpp_Species] cpp_reactants
#     for sp in reactants:
#         cpp_reactants.push_back(deref((<Species> sp).thisptr))
#     return context.count_rrmatches(deref(pttrn.thisptr), cpp_reactants)
# 
# def rrgenerate(ReactionRule pttrn, reactants):
#     """rrgenerate(pttrn, reactants) -> [Species]
# 
#     Generate a list of products from the given list of reactants.
# 
#     Parameters
#     ----------
#     pttrn : ReactionRule
#         A pattern.
#     reactants : list
#         A list of ``Species``. The order of ``reactants`` is respected.
# 
#     Returns
#     -------
#     list:
#         A list of products. The size of the list is equal to the number of matches.
#         Each element of the list is a list of ``Species``.
# 
#     Notes
#     -----
#     Rather use ``ReactionRule.generate``.
# 
#     """
#     cdef vector[Cpp_Species] cpp_reactants
#     for sp in reactants:
#         cpp_reactants.push_back(deref((<Species> sp).thisptr))
#     cdef vector[vector[Cpp_Species]] cpp_products_list = \
#         context.rrgenerate(deref(pttrn.thisptr), cpp_reactants)
#     cdef vector[vector[Cpp_Species]].iterator it1 = cpp_products_list.begin()
#     cdef vector[Cpp_Species].iterator it2
#     retval = []
#     while it1 != cpp_products_list.end():
#         retval.append([])
#         it2 = deref(it1).begin()
#         while it2 != deref(it1).end():
#             retval[-1].append(Species_from_Cpp_Species(address(deref(it2))))
#             inc(it2)
#         inc(it1)
#     return retval
