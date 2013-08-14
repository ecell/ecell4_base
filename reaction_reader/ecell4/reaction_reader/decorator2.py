import types
import numbers
import species
import copy
import decorator
import parseobj
import functools


def generate_Species(obj):
    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        sp = species.Species()
        for elem in obj._elements():
            su = species.Subunit(elem.name)
            if elem.args is not None:
                for mod in elem.args:
                    if (not (isinstance(mod, parseobj.ParseObj)
                            or isinstance(mod, parseobj.AnyCallable))
                        or mod._size() != 1):
                        raise RuntimeError, (
                            "invalid argument [%s] found." % (mod))
                    arg = mod._elements()[0]
                    name, binding = arg.name, arg.modification
                    if binding is None:
                        su.add_modification(name, "", "")
                    else:
                        binding = str(binding)
                        if not (binding.isdigit() or binding == ""
                            or binding[0] == "_"):
                            raise RuntimeError, (
                                "invalid binding [%s] given." % (binding))
                        su.add_modification(name, "", binding)
            if elem.kwargs is not None:
                for name, value in elem.kwargs.items():
                    if (not (isinstance(value, parseobj.ParseObj)
                            or isinstance(value, parseobj.AnyCallable))
                        or value._size() != 1):
                        raise RuntimeError, (
                            "invalid argument [%s] found." % (value))
                    arg = value._elements()[0]
                    state, binding = str(arg.name), arg.modification
                    if binding is None:
                        su.add_modification(name, state, "")
                    else:
                        binding = str(binding)
                        if not (binding.isdigit() or binding == ""
                            or binding[0] == "_"):
                            raise RuntimeError, (
                                "invalid binding [%s] given." % (binding))
                        su.add_modification(name, state, binding)
            sp.add_subunit(su)
        return (sp, )
    elif isinstance(obj, parseobj.InvExp):
        return (None, )
    elif isinstance(obj, parseobj.AddExp):
        subobjs = obj._elements()
        return tuple(generate_Species(subobj)[0] for subobj in subobjs)

    raise RuntimeError, 'invalid expression; "%s" given' % str(obj)

def generate_ReactionRule(lhs, rhs, opts):
    # if len(lhs) == 0:
    #     if len(rhs) != 1:
    #         raise RuntimeError, (
    #             "the number of products must be 1; %d given" % len(rhs))
    #     return ecell4.core.create_synthesis_reaction_rule(rhs[0], k)
    if len(lhs) == 0 or len(lhs) == 1 or len(lhs) == 2:
        return species.ReactionRule(lhs, rhs, opts)
    raise RuntimeError, (
        "the number of reactants must be less than 3; %d given" % len(lhs))

def generate_Option(opt):
    # if not (isinstance(opt, parseobj.AnyCallable)
    #     or isinstance(opt, parseobj.ParseObj)):
    #     raise RuntimeError

    if opt._size() != 1:
        raise RuntimeError

    elem = opt._elements()[0]
    if elem.name == "IncludeReactants" or elem.name == "ExcludeReactants":
        if not (len(elem.args) == 2
            and type(elem.args[0]) == int
            and (isinstance(elem.args[1], parseobj.AnyCallable)
                or isinstance(elem.args[1], parseobj.ParseObj))):
            raise RuntimeError

        if isinstance(elem.args[1], parseobj.ParseObj):
            raise RuntimeError, "only a subunit name is allowed here."

        pttrn = elem.args[1]._elements()[0].name
        if elem.name == "ExcludeReactants":
            return (species.ExcludeReactants(elem.args[0], pttrn),
                species.ExcludeProducts(elem.args[0], pttrn))
        elif elem.name == "IncludeReactants":
            return (species.IncludeReactants(elem.args[0], pttrn),
                species.IncludeProducts(elem.args[0], pttrn))
    elif elem.name == "IncludeProducts" or elem.name == "ExcludeProducts":
        if not (len(elem.args) == 2
            and type(elem.args[0]) == int
            and (isinstance(elem.args[1], parseobj.AnyCallable)
                or isinstance(elem.args[1], parseobj.ParseObj))):
            raise RuntimeError

        if isinstance(elem.args[1], parseobj.ParseObj):
            raise RuntimeError, "only a subunit name is allowed here."

        pttrn = elem.args[1]._elements()[0].name
        if elem.name == "ExcludeProducts":
            return (species.ExcludeProducts(elem.args[0], pttrn),
                species.ExcludeReactants(elem.args[0], pttrn))
        elif elem.name == "IncludeProducts":
            return (species.IncludeProducts(elem.args[0], pttrn),
                species.IncludeReactants(elem.args[0], pttrn))
    else:
        # raise RuntimeError
        return (opt, None)
    return (opt, opt)

def generate_Options1(opts):
    retval = []
    for opt in opts:
        if (isinstance(opt, parseobj.AnyCallable)
            or isinstance(opt, parseobj.ParseObj)):
            lhs, rhs = generate_Option(opt)
            if lhs is not None:
                retval.append(lhs)
        # elif isinstance(opt, numbers.Number):
        #     retval.append(opt)
        # else:
        #     raise RuntimeError, "an invalid option [%s] given." % (opt)
        retval.append(opt)
    return retval

def generate_Options2(opts):
    retval1, retval2 = [], []
    for opt in opts:
        if (isinstance(opt, parseobj.AnyCallable)
            or isinstance(opt, parseobj.ParseObj)):
            lhs, rhs = generate_Option(opt)
            if lhs is not None:
                retval1.append(lhs)
            if rhs is not None:
                retval2.append(rhs)
        elif ((isinstance(opt, types.ListType)
            or isinstance(opt, types.TupleType))
            and len(opt) == 2):
            # if (isinstance(opt[0], numbers.Number)
            #     and isinstance(opt[1], numbers.Number)):
            #     raise RuntimeError
            retval1.append(opt[0])
            retval2.append(opt[1])
        else:
            raise RuntimeError, "an invalid option [%s] given." % (opt)
    return retval1, retval2

class SpeciesAttributesCallback(decorator.Callback):

    def __init__(self, *args):
        decorator.Callback.__init__(self)

        self.bitwise_operations = []

    def get(self):
        return copy.copy(self.bitwise_operations)

    def notify_bitwise_operations(self, obj):
        if not isinstance(obj, parseobj.OrExp):
            raise RuntimeError, 'an invalid object was given [%s]' % (repr(obj))
        elif len(obj._elements()) != 2:
            raise RuntimeError, 'only one attribute is allowed. [%d] given' % (
                len(obj._elements()))

        lhs, rhs = obj._elements()

        species_list = generate_Species(lhs)
        if len(species_list) != 1:
            raise RuntimeError, (
                'only a single species must be given; %d given'
                % len(species_list))

        sp = species_list[0]
        if sp is None:
            raise RuntimeError, "no species given [%s]" % (repr(obj))

        self.bitwise_operations.append((sp, rhs))

    def notify_comparisons(self, obj):
        raise RuntimeError, (
            'ReactionRule definitions are not allowed'
            + ' in "species_attributes"')

class ReactionRulesCallback(decorator.Callback):

    def __init__(self):
        decorator.Callback.__init__(self)

        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_comparisons(self, obj):
        if not isinstance(obj, parseobj.CmpExp):
            raise RuntimeError, 'an invalid object was given [%s]' % (repr(obj))
        elif isinstance(obj, parseobj.NeExp):
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)

        lhs, rhs = obj._lhs, obj._rhs

        if isinstance(lhs, parseobj.OrExp):
            lhs = lhs._elements()[0]

        if isinstance(rhs, parseobj.OrExp):
            opts = rhs._elements()[1: ]
            rhs = rhs._elements()[0]
        else:
            opts = []

        lhs, rhs = generate_Species(lhs), generate_Species(rhs)
        lhs = tuple(sp for sp in lhs if sp is not None)
        rhs = tuple(sp for sp in rhs if sp is not None)

        if isinstance(obj, parseobj.EqExp) or isinstance(obj, parseobj.NeExp):
            opts1, opts2 = generate_Options2(opts)
            self.comparisons.append(generate_ReactionRule(lhs, rhs, opts1))
            self.comparisons.append(generate_ReactionRule(rhs, lhs, opts2))
        elif isinstance(obj, parseobj.GtExp):
            opts = generate_Options1(opts)
            self.comparisons.append(generate_ReactionRule(lhs, rhs, opts))
        else:
            raise RuntimeError, 'an invalid object was given [%s]' % (repr(obj))

species_attributes = functools.partial(decorator.parse_decorator, SpeciesAttributesCallback)
reaction_rules = functools.partial(decorator.parse_decorator, ReactionRulesCallback)
