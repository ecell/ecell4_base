import copy
import types
import numbers
import warnings
import functools

import parseobj

import ecell4.core


#============================================================
# Species Meta Infomation CLASS XXX
class Subunit:
    def __init__(self, name = None):
        self.name = name
        self.modification_list = []
        self.species = None
    def add_modification(self, new_modification):
        self.modification_list.append(new_modification)
    def get_modification(self):
        return self.modification_list

    # XXX Following 2 functions( (set|get)_species ) may not be necessary.
    #   Because it makes cross refenrence between Subunit and Meta_Species, 
    #       the structure of class dependency will be complicated.
    def set_species(self, sp):
        self.species = sp
    def get_species(self):
        return self.species
    def enum_binding_subunit(self):
        ret = []
        for mod in self.modification_list:
            for binding in mod.get_binding():
                ret.append(binding.subunit)
        return ret
    def name(self):
        pass

class Modification:
    def __init__(self, subunit, name, attribute = None):
        self.subunit = subunit  # reference 
        self.name = name
        self.binding = []   # references to binding partner's Modification object.
        self.attribute = attribute 

    def set_binding(self, substrates):
        self.binding = [s for s in substrates if s is not self]
    def get_binding(self):
        return self.binding
    def get_subunit(self):
        return self.subunit
    def name(self):
        return self.name

class Meta_Species(ecell4.core.Species):
    def __init__(self, name):
        self.subunit_list = []
        ecell4.core.Species.__init__(self, name)
    def add_subunit(self, sub):
        self.subunit_list.append(sub)
        sub.set_species(self)
    def get_subunit(self):
        return self.subunit_list

    def dot_output(self):
        acc = []
        print "digraph %s {" % self.name()
        print "\tgraph [label = \"%s\", labelloc = t];" % self.name()
        for sub in self.subunit_list:
            for s in sub.enum_binding_subunit():
                print "\t\"%s\" -> \"%s\";" % (sub.name, s.name)
        print "}"


#============================================================

def generate_Species2(obj):
    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        # COMPLEX
        elems = obj._get_elements()
        msp = Meta_Species(elems[0].name)
        correct_binding_dict = {}

        for e in elems:
            s = Subunit(e)
            msp.add_subunit(s)
            if not isinstance(e.args, types.NoneType):
                for arg in e.args:
                    parse_elem_list = arg._get_elements()
                    for parse_elem in parse_elem_list:
                        modification_info = Modification(s, parse_elem.name)
                        s.add_modification(modification_info)
                        if parse_elem.modification in correct_binding_dict:
                            correct_binding_dict[ parse_elem.modification ].append(modification_info)
                        else:
                            correct_binding_dict[ parse_elem.modification ] = [ modification_info ]

            if not isinstance(e.args, types.NoneType):
                for k, v in e.kwargs.iteritems():
                    parse_elem_list = v._get_elements()
                    for parse_elem in parse_elem_list:
                        modification_info = Modification(s, parse_elem.name, k)
                        s.add_modification(modification_info)
                        if parse_elem.modification in correct_binding_dict:
                            correct_binding_dict[ parse_elem.modification ].append(modification_info)
                        else:
                            correct_binding_dict[ parse_elem.modification ] = [ modification_info ]
        for index, array in correct_binding_dict.iteritems():
            for modification_iter in array:
                modification_iter.set_binding(array)
        msp.dot_output()
        return [msp]
        
    elif isinstance(obj, parseobj.ParseObjSet):
        subobjs = obj._get_objects()
        species_list = []
        for s in subobjs:
            species_list.extend( generate_Species2(s) )
        return species_list
        

def generate_Species(obj):
    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        elems = obj._get_elements()
        if len(elems) != 1:
            raise NotImplementedError, (
                'complex is not allowed yet; "%s"' % str(obj))
        if (elems[0].args is not None
            or elems[0].kwargs is not None
            or elems[0].key is not None
            or elems[0].modification is not None):
            raise NotImplementedError, (
                'modification is not allowed yet; "%s"' % str(obj))
        if elems[0].inv:
            return ((None, elems[0].param), )
        else:
            return ((ecell4.core.Species(elems[0].name), elems[0].param), )
    elif isinstance(obj, parseobj.ParseObjSet):
        subobjs = obj._get_objects()
        return tuple(generate_Species(subobj)[0] for subobj in subobjs)
    raise RuntimeError, 'invalid expression; "%s" given' % str(obj)


def generate_ReactionRule2(lhs, rhs, k=0.0):
    #arguments (: lhs rhs) are lists of Meta_Species objects.
    if len(lhs) == 0:
        if len(rhs) != 1:
            raise RuntimeError, (
                "the number of products must be 1; %d given" % len(rhs))
        return ecell4.core.create_synthesis_reaction_rule(rhs[0], k)
    elif len(lhs) == 1:
        if len(rhs) == 0:
            return ecell4.core.create_degradation_reaction_rule(lhs[0], k)
        elif len(rhs) == 1:
            return ecell4.core.create_unimolecular_reaction_rule(
                lhs[0], rhs[0], k)
        elif len(rhs) == 2:
            return ecell4.core.create_unbinding_reaction_rule(
                lhs[0], 
                rhs[0], 
                rhs[1], k)
        else:
            raise RuntimeError, (
                "the number of products must be less than 3; %d given"
                % len(rhs))
    elif len(lhs) == 2:
        if len(rhs) == 1:
            return ecell4.core.create_binding_reaction_rule(
                lhs[0], 
                lhs[1], 
                rhs[0], k)
        else:
            raise RuntimeError, (
                "the number of products must be 1; %d given" % len(rhs))
    raise RuntimeError, (
        "the number of reactants must be less than 3; %d given" % len(lhs))


def generate_ReactionRule(lhs, rhs, k=0.0):
    if len(lhs) == 0:
        if len(rhs) != 1:
            raise RuntimeError, (
                "the number of products must be 1; %d given" % len(rhs))
        return ecell4.core.create_synthesis_reaction_rule(rhs[0], k)
    elif len(lhs) == 1:
        if len(rhs) == 0:
            return ecell4.core.create_degradation_reaction_rule(lhs[0], k)
        elif len(rhs) == 1:
            return ecell4.core.create_unimolecular_reaction_rule(
                lhs[0], rhs[0], k)
        elif len(rhs) == 2:
            return ecell4.core.create_unbinding_reaction_rule(
                lhs[0], rhs[0], rhs[1], k)
        else:
            raise RuntimeError, (
                "the number of products must be less than 3; %d given"
                % len(rhs))
    elif len(lhs) == 2:
        if len(rhs) == 1:
            return ecell4.core.create_binding_reaction_rule(
                lhs[0], lhs[1], rhs[0], k)
        else:
            raise RuntimeError, (
                "the number of products must be 1; %d given" % len(rhs))
    raise RuntimeError, (
        "the number of reactants must be less than 3; %d given" % len(lhs))

class SpeciesAttributesCallback(object):

    def __init__(self, *args):
        self.keys = None
        if len(args) > 0:
            for key in args:
                if not isinstance(key, types.StringType):
                    raise RuntimeError, 'non string key "%s" was given' % key
            self.keys = args

        self.bitwise_operations = []

    def get(self):
        return copy.copy(self.bitwise_operations)

    def notify_unary_operations(self, optr, target):
        pass

    def notify_bitwise_operations(self, optr, lhs, rhs):
        if optr != '|':
            raise RuntimeError, 'operator "%s" not allowed' % optr

        species_list = generate_Species(lhs)
        if len(species_list) != 1:
            raise RuntimeError, (
                'only a single species must be given; %d given'
                % len(species_list))

        sp, _ = species_list[0]
        if sp is None:
            raise RuntimeError, 'never use "~" in "species_attributes"'

        if self.keys is None:
            if not isinstance(rhs, types.DictType):
                raise RuntimeError, (
                    'parameter must be given as a dict; "%s" given'
                    % str(rhs))
            for key, value in rhs.items():
                if not (isinstance(key, types.StringType)
                    and isinstance(value, types.StringType)):
                    raise RuntimeError, (
                        'attributes must be given as a pair of strings;'
                        + ' "%s" and "%s" given'
                        % (str(key), str(value)))
                sp.set_attribute(key, value)
        else:
            if not (isinstance(rhs, types.TupleType)
                and isinstance(rhs, types.ListType)):
                if len(self.keys) == 1:
                    rhs = (rhs, )
                else:
                    raise RuntimeError, (
                        'parameters must be given as a tuple or list; "%s" given'
                        % str(rhs))
            if len(rhs) != len(self.keys):
                raise RuntimeError, (
                    'the number of parameters must be %d; %d given'
                    % (len(self.keys), len(rhs)))
            else:
                for key, value in zip(self.keys, rhs):
                    if not isinstance(value, types.StringType):
                        raise RuntimeError, (
                            'paramter must be given as a string; "%s" given'
                            % str(value))
                    sp.set_attribute(key, value)

        self.bitwise_operations.append(sp)

    def notify_comparisons(self, optr, lhs, rhs):
        raise RuntimeError, (
            'ReactionRule definitions are not allowed'
            + ' in "species_attributes"')

class ReactionRulesCallback(object):

    def __init__(self):
        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_unary_operations(self, optr, target):
        pass

    def notify_bitwise_operations(self, optr, lhs, rhs):
        pass

    def notify_comparisons(self, optr, lhs, rhs):
        if optr == "!=":
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)

        lhs, rhs = generate_Species(lhs), generate_Species(rhs)
        lhs = tuple(sp for sp, param in lhs if sp is not None)
        k = rhs[-1][1]
        rhs = tuple(sp for sp, param in rhs if sp is not None)

        if optr == "!=" or optr == "==":
            if not (isinstance(k, types.ListType)
                or isinstance(k, types.TupleType)):
                raise RuntimeError, (
                    'parameter must be a list or tuple with length 2; "%s" given'
                    % str(k))
            elif len(k) != 2:
                raise RuntimeError, (
                    "parameter must be a list or tuple with length 2;"
                    + " length %d given" % len(k))
            elif not (isinstance(k[0], numbers.Number)
                and isinstance(k[1], numbers.Number)):
                raise RuntimeError, (
                    'parameters must be given as a list or tuple of numbers;'
                    + ' "%s" given' % str(k))
            self.comparisons.append(generate_ReactionRule(lhs, rhs, k[0]))
            self.comparisons.append(generate_ReactionRule(rhs, lhs, k[1]))
        elif optr == ">":
            if k is None:
                raise RuntimeError, 'no parameter is specified'
            elif not isinstance(k, numbers.Number):
                raise RuntimeError, (
                    'parameter must be given as a number; "%s" given' % str(k))
            self.comparisons.append(generate_ReactionRule(lhs, rhs, k))
        else:
            raise RuntimeError, 'operator "%s" not allowed' % optr


class JustParseCallback(object):

    def __init__(self):
        self.comparisons = []
        self.kinetic_parameters = []
        # self.bitwise_optrs = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_unary_operations(self, optr, target):
        pass

    def notify_bitwise_operations(self, optr, lhs, rhs):
        if isinstance(rhs, types.TupleType):
            self.kinetic_parameters = list(rhs)# kon, koff
        else:
            self.kinetic_parameters = [rhs]

    def notify_comparisons(self, optr, lhs, rhs):
        if optr == "!=":
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)

        # After calling next sentence(2 times of generate_Species2) , 
        #   the left side variables (lhs and rhs)  will be the array of Meta_Species.
        lhs, rhs = generate_Species2(lhs), generate_Species2(rhs)

        if optr == "==" or optr == "!=":
            if len(self.kinetic_parameters) != 2:
                raise RuntimeError("The number of kinetic parameters is invalid.")
            # reversible reaction
            self.comparisons.append(generate_ReactionRule2(lhs, rhs, self.kinetic_parameters[0]))
            self.comparisons.append(generate_ReactionRule2(rhs, lhs, self.kinetic_parameters[1]))
        elif optr == ">":
            # irreversible reaction
            self.comparisons.append(generate_ReactionRule2(lhs, rhs, self.kinetic_parameters[0]))
            pass
        else:
            raise RuntimeError, 'operator "%s" not allowed' % optr

class Callback(object):
    """callback before the operations"""

    def __init__(self):
        # self.bitwise_operations = []
        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_unary_operations(self, optr, target):
        pass

    def notify_bitwise_operations(self, optr, lhs, rhs):
        # self.bitwise_operations.append((optr, lhs, rhs))
        pass

    def notify_comparisons(self, optr, lhs, rhs):
        if optr == "!=":
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)
        self.comparisons.append((optr, lhs, rhs))

def parse_decorator(callback_class, func):
    @functools.wraps(func)
    def wrapped(*args, **kwargs):
        cache = callback_class()
        vardict = copy.copy(func.func_globals)
        #print vardict.keys()
        for k in func.func_code.co_names:
            if (not k in vardict.keys()
                and not k in dir(vardict['__builtins__'])): # is this enough?
                vardict[k] = parseobj.AnyCallable(cache, k)
        g = types.FunctionType(func.func_code, vardict)
        with warnings.catch_warnings():
            # warnings.simplefilter("always")
            g(*args, **kwargs)
        return cache.get()
    return wrapped

# reaction_rules = functools.partial(parse_decorator, Callback)
reaction_rules = functools.partial(parse_decorator, ReactionRulesCallback)
species_attributes = functools.partial(parse_decorator, SpeciesAttributesCallback)
just_parse = functools.partial(parse_decorator, JustParseCallback)

def species_attributes_with_keys(*args):
    def create_callback():
        return SpeciesAttributesCallback(*args)
    return functools.partial(parse_decorator, create_callback)
