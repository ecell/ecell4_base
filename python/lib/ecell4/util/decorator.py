import copy
import numbers
import functools

from . import parseobj
from .decorator_base import Callback, JustParseCallback, ParseDecorator

from .model_parser import (
    generate_species,
    generate_list_of_species_with_coefficient,
    generate_reaction_rule_options,
    generate_reaction_rule,
    as_quantity
    )

import ecell4.core

SPECIES_ATTRIBUTES = []
REACTION_RULES = []

class SpeciesAttributesCallback(Callback):

    def __init__(self):
        Callback.__init__(self)

        self.bitwise_operations = []

    def get(self):
        return copy.copy(self.bitwise_operations)

    def set(self):
        global SPECIES_ATTRIBUTES
        SPECIES_ATTRIBUTES.extend(self.bitwise_operations)

    def notify_bitwise_operations(self, obj):
        self.bitwise_operations.extend(self.generate(obj))

    @staticmethod
    def generate(obj):
        # attribute_dimensionality = {'D': '[length]**2/[time]', 'radius': '[length]'}

        if not isinstance(obj, parseobj.OrExp):
            raise TypeError('An invalid object was given [{}]'.format(repr(obj)))

        elems = obj._elements()
        rhs = elems[-1]
        if isinstance(rhs, parseobj.ExpBase):
            return ()
        elif not isinstance(rhs, dict):
            raise TypeError('parameter must be given as a dict; "{}" given'.format(type(rhs)))

        species_list = []
        for lhs in elems[: -1]:
            sp = generate_species(lhs)
            for key, value in rhs.items():
                if not isinstance(key, str):
                    raise TypeError(
                        "Attribute key must be string."
                        " '{}' was given [{}].".format(type(key).__name__, key))

                value = as_quantity(value)
                if not isinstance(value, (numbers.Real, str, bool, ecell4.core.Quantity)):
                    raise TypeError(
                        "Attribute value must be int, float, string, boolean or Quantity."
                        " '{}' was given [{}].".format(type(value).__name__, value))
                sp.set_attribute(key, value)
            species_list.append(sp)
        return species_list

    def notify_comparisons(self, obj):
        raise RuntimeError(
            'ReactionRule definitions are not allowed'
            + ' in "species_attributes"')

class ReactionRulesCallback(Callback):

    def __init__(self):
        Callback.__init__(self)

        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def set(self):
        global REACTION_RULES
        REACTION_RULES.extend(self.comparisons)

    def notify_comparisons(self, obj):
        self.comparisons.extend(self.generate(obj))

    @staticmethod
    def generate(obj):
        if not isinstance(obj, (parseobj.EqExp, parseobj.GtExp, parseobj.LtExp)):
            raise RuntimeError('An ivalid operation was given [{!r}].'.format(obj))

        lhs, rhs = obj._lhs, obj._rhs

        if isinstance(obj._lhs, parseobj.OrExp):
            lhs = obj._lhs._elements()[0]
        else:
            lhs = obj._lhs

        if not isinstance(obj._rhs, parseobj.OrExp):
            raise RuntimeError(
                "A right-hand-side is ill-formed. OrExp must be given [{!r}].".format(obj._rhs))
        elif len(obj._rhs._elements()) == 0:
            raise RuntimeError('No product was given in the right-hand-side.')
        else:
            rhs = obj._rhs._elements()[0]
            opts = obj._rhs._elements()[1: ]

        reactants = generate_list_of_species_with_coefficient(lhs)
        products = generate_list_of_species_with_coefficient(rhs)
        reactants = tuple((sp, coef) for sp, coef in reactants if coef != 0)
        products = tuple((sp, coef) for sp, coef in products if coef != 0)

        opts = generate_reaction_rule_options(opts)
        if 'k' not in opts.keys():
            raise RuntimeError('No kinetic rate or law was given.')
        params = opts['k']

        if isinstance(obj, parseobj.EqExp):
            if not isinstance(params, (tuple, list)):
                raise RuntimeError(
                    "Parameter must be list or tuple."
                    " '{}' was given [{}].".format(type(params).__name__, params))
            elif len(params) != 2:
                raise RuntimeError(
                    "Parameter must have size, 2."
                    " '{}' was given [{}].".format(len(params), params))

            return (generate_reaction_rule(reactants, products, as_quantity(params[0]), opts.get('policy')),
                    generate_reaction_rule(products, reactants, as_quantity(params[1]), opts.get('policy')))
        elif isinstance(obj, parseobj.GtExp):
            return (generate_reaction_rule(reactants, products, as_quantity(params), opts.get('policy')), )
        elif isinstance(obj, parseobj.LtExp):
            return (generate_reaction_rule(products, reactants, as_quantity(params), opts.get('policy')), )

def get_model(is_netfree=False, without_reset=False, seeds=None, effective=False):
    """
    Generate a model with parameters in the global scope, ``SPECIES_ATTRIBUTES``
    and ``REACTIONRULES``.

    Parameters
    ----------
    is_netfree : bool, optional
        Return ``NetfreeModel`` if True, and ``NetworkModel`` if else.
        Default is False.
    without_reset : bool, optional
        Do not reset the global variables after the generation if True.
        Default is False.
    seeds : list, optional
        A list of seed ``Species`` for expanding the model.
        If this is not None, generate a ``NetfreeModel`` once, and return a
        ``NetworkModel``, which is an expanded form of that with the given seeds.
        Default is None.
    effective : bool, optional
        See ``NetfreeModel.effective`` and ``Netfree.set_effective``.
        Only meaningfull with option ``is_netfree=True``.
        Default is False

    Returns
    -------
    model : NetworkModel, NetfreeModel

    """
    try:
        if seeds is not None or is_netfree:
            m = ecell4.core.NetfreeModel()
        else:
            m = ecell4.core.NetworkModel()

        for sp in SPECIES_ATTRIBUTES:
            m.add_species_attribute(sp)
        for rr in REACTION_RULES:
            m.add_reaction_rule(rr)

        if not without_reset:
            reset_model()

        if seeds is not None:
            return m.expand(seeds)

        if isinstance(m, ecell4.core.NetfreeModel):
            m.set_effective(effective)
    except Exception as e:
        reset_model()
        raise e

    return m

def reset_model():
    """
    Reset all values, ``SPECIES_ATTRIBUTES`` and ``REACTIONRULES``,
    in the global scope.

    """
    global SPECIES_ATTRIBUTES
    global REACTION_RULES

    SPECIES_ATTRIBUTES = []
    REACTION_RULES = []

reaction_rules = functools.partial(ParseDecorator, ReactionRulesCallback)
species_attributes = functools.partial(ParseDecorator, SpeciesAttributesCallback)
