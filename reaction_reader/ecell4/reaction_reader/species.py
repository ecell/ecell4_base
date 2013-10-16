import copy
import itertools
import sys

import options


label_subunit = lambda x: "subunit%s" % x
label_binding = lambda x: "binding%s" % x
label_domain = lambda x, y: "_%s__%s" % (x, y)

def is_label_domain(key, subunit):
    header = "_%s__" % subunit
    return (len(key) > len(header) and key[: len(header)] == header)

class Species(object):

    def __init__(self):
        self.subunits = []

        self.conditions = None

    def num_bindings(self):
        labels, n = [], 0
        for subunit in self.subunits:
            for mod, (state, binding) in subunit.modifications.items():
                if binding == "":
                    continue
                elif binding[0] == "_":
                    if len(binding) != 1:
                        raise RuntimeError, "[%s] not supported yet." % binding
                    n += 1
                else:
                    labels.append(int(binding))
        return n + len(set(labels))

    def get_binding_stride(self):
        retval = 0
        for subunit in self.subunits:
            for mod, (state, binding) in subunit.modifications.items():
                if binding != "" and binding[0] != "_":
                    retval = max(retval, int(binding))
        return retval

    def update_indices(self):
        for i, subunit in enumerate(self.subunits):
            subunit.index = i

    def add_subunit(self, subunit):
        subunit.index = len(self.subunits)
        self.subunits.append(subunit)

    def get_subunit_list(self):
        return self.subunits

    def num_subunits(self, pttrn=None):
        if pttrn is None:
            return len(self.subunits)

        retval = 0
        for su in self.subunits:
            if su.name == pttrn:
                retval += 1
        return retval

    def generate_conditions(self, stride=0):
        conditions = []
        for i, subunit in enumerate(self.subunits):
            conditions.extend(
                subunit.generate_conditions(label_subunit(stride + i)))
        conditions.sort(key=lambda x: x.priority)
        return conditions

    def match(self, sp, contexts=None, stride=None):
        if contexts is None:
            contexts = Contexts()
            contexts.initialize()
        elif len(contexts) == 0:
            return contexts

        if stride is None:
            stride = contexts.num_subunits()

        if self.conditions is None:
            self.conditions = self.generate_conditions(stride)

        for condition in self.conditions:
            contexts = condition.match(sp, contexts)
            if len(contexts) == 0:
                break

        # conditions = self.generate_conditions(stride)

        # for condition in conditions:
        #     contexts = condition.match(sp, contexts)
        #     if len(contexts) == 0:
        #         break

        contexts.clear_locals()
        return contexts

    def sort(self):
        cmpsu = CmpSubunit(self)
        cmpsu.sort()

    def __str__(self):
        # self.sort() #XXX: check if it's already sorted or not
        return ".".join([str(subunit) for subunit in self.subunits])

    def __repr__(self):
        return '<"%s">' % (str(self))

    def __eq__(self, rhs):
        if len(self.subunits) != len(rhs.subunits):
            # quick filtering
            return False
        else:
            return (len(self.match(rhs)) > 0 and len(rhs.match(self)) > 0)

class Subunit(object):

    def __init__(self, name):
        self.name = name
        self.modifications = {}
        self.exclusions = []
        self.domain_classes = {}
        self.commutatives = Commutatives()

        self.index = None #XXX

    def get_name(self):
        return self.name

    def generate_conditions(self, key):
        conditions = []
        conditions.append(SubunitContainingCondition(key, self.name))
        for mod in self.exclusions:
            conditions.append(
                ExcludedModificationCondition(key, mod))

        for mod, value in self.domain_classes.items():
            conditions.append(
                DomainClassCondition(key, mod, value))

        for value in self.commutatives.as_sets():
            conditions.append(
                CommutativeCondition(key, value))

        for mod, (state, binding) in self.modifications.items():
            conditions.append(
                ModificationNameCondition(key, mod))
            conditions.append(
                ModificationBindingCondition(key, mod, binding))
            if state != "":
                conditions.append(
                    ModificationStateCondition(key, mod, state))
        return conditions

    def add_modification(self, mod, state="", binding=""):
        if isinstance(binding, str) and len(binding) > 1 and binding[0] == "_":
            raise RuntimeError, "A binding label [%s] not allowed" % (binding)

        self.modifications[mod] = (state, str(binding))

    def get_modifications_list(self):
        return self.modifications

    def add_exclusion(self, mod):
        if not mod in self.exclusions:
            self.exclusions.append(mod)

    def add_domain_class(self, mod, value):
        if (not (isinstance(value, list) or isinstance(value, tuple))
            or len(value) == 0):
            raise ValueError, "Invalid argument [%s] given." % str(value)

        self.domain_classes[mod] = value

    def set_commutative(self, mods):
        self.commutatives.set_commutative(*mods)

    def __str__(self):
        mods1 = ["~%s" % (mod) for mod in self.exclusions]

        mods2, mods3 = [], []
        for mod, (state, binding) in self.modifications.items():
            if state == "":
                if binding != "":
                    mods2.append("%s^%s" % (mod, binding))
                else:
                    mods2.append(mod)
            elif binding == "":
                mods3.append("%s=%s" % (mod, state))
            else:
                mods3.append("%s=%s^%s" % (mod, state, binding))

        mods4 = ["%s=[%s]" % (mod, ",".join([str(elem) for elem in value]))
                for mod, value in self.domain_classes.items()]

        mods5 = ["(%s)" % (",".join(subset))
            for subset in self.commutatives.as_sets()]

        mods1.sort()
        mods2.sort()
        mods3.sort()
        mods4.sort()
        mods5.sort()

        labels = ",".join(itertools.chain(mods1, mods2, mods5, mods3, mods4))
        if labels == "":
            return self.name
        else:
            return "%s(%s)" % (self.name, labels)

    def __repr__(self):
        return '<"%s">' % (str(self))

class Commutatives(object):

    def __init__(self):
        self.__indices = {}

    def __len__(self):
        return len(self.__indices)

    def keys(self):
        return self.__indices.keys()

    def get_commutatives(self, key1):
        idx1 = self.__indices.get(key1)
        if idx1 is None:
            return []
        else:
            return [key2 for key2, idx2 in self.__indices.items()
                if idx2 == idx1]

    def set_commutative(self, *keys):
        if len(keys) == 0:
            return

        collection, newset = [], []
        for key in keys:
            if key in self.__indices.keys():
                collection.append(self.__indices[key])
            else:
                newset.append(key)

        if len(collection) == 0:
            newidx = (max(self.__indices.values()) + 1
                if len(self.__indices) > 0 else 1)
        else:
            newidx = min(collection)
            for key in self.__indices.keys():
                if self.__indices[key] in collection:
                    self.__indices[key] = newidx

        for key in newset:
            self.__indices[key] = newidx

    def is_commutative(self, *keys):
        if len(keys) < 2:
            raise RuntimeError, "at least two keys must be given"

        idx1 = self.__indices.get(keys[0])
        if idx1 is None:
            # raise RuntimeError, "an invalid key [%s] given" % (keys[0])
            return False

        for key in keys[1: ]:
            idx2 = self.__indices.get(key)
            if idx2 is None:
                # raise RuntimeError, "an invalid key [%s] given" % (key)
                return False
            elif idx2 != idx1:
                return False
        return True

    def is_comprised(self, com):
        if len(com) < len(self):
            return False

        for subset in self.as_sets():
            if not com.is_commutative(*subset):
                return False
        return True

    def as_sets(self):
        return commutative_generator(self.__indices)

def commutative_generator(indices):
    i, done, keys = 0, [], indices.keys()

    for i in range(len(keys)):
        value = indices[keys[i]]
        if value in done:
            continue
        done.append(value)
        retval = [key for key in keys[i: ] if indices[key] == value]
        retval.sort()
        yield tuple(retval)

def check_connectivity(src, markers=[]):
    adjacencies = {}
    tmp = {}
    for i, su in enumerate(src.subunits):
        adjacencies[i] = []
        for mod, (state, binding) in su.modifications.items():
            if binding == "":
                continue

            if binding in tmp.keys():
                if tmp[binding] is None:
                    raise RuntimeError, "[%s] duplicated in [%s:%d]" % (
                        binding, src, src.get_binding_stride())

                adjacencies[i].append(tmp[binding])
                adjacencies[tmp[binding]].append(i)
                tmp[binding] = None
            else:
                tmp[binding] = i
    for binding, value in tmp.items():
        if value is not None:
            raise RuntimeError, "no target for [%s] in [%s]" % (binding, src)

    L = list(range(len(src.subunits)))
    Ks = []
    while len(L) != 0:
        K = [L.pop()]
        reconnect(L, K, adjacencies)
        Ks.append(K)

    if len(Ks) == 0:
        return (None, [None for _ in markers])
    elif len(Ks) == 1:
        return ((src, ), [0 for _ in markers])
    else:
        products, correspondence = [], [None for _ in markers]
        for K in Ks:
            sp = Species()
            for i in K:
                sp.add_subunit(copy.deepcopy(src.subunits[i]))
                if src.subunits[i].index in markers:
                    correspondence[markers.index(src.subunits[i].index)] = len(products)
            sp.update_indices()
            products.append(sp)
        return (tuple(products), correspondence)

def reconnect(L, K, adjacencies):
    src = K[-1]
    for i in adjacencies[src]:
        if not i in K:
            if i in L:
                L.remove(i)
            K.append(i)
            reconnect(L, K, adjacencies)

def concatenate_species(*species_list):
    retval, stride = Species(), 0
    for sp in species_list:
        for su in sp.subunits:
            newsu = Subunit(su.name)
            for mod in su.exclusions:
                newsu.add_exclusion(mod)
            for key, value in su.domain_classes.items():
                newsu.add_domain_class(key, value)
            for subset in su.commutatives.as_sets():
                newsu.set_commutative(subset)
            for mod, (state, binding) in su.modifications.items():
                if binding != "" and binding[0] != "_":
                    if not binding.isdigit():
                        raise RuntimeError
                    binding = int(binding) + stride
                newsu.add_modification(mod, state, binding)
            retval.add_subunit(newsu)
        stride += sp.get_binding_stride()
    retval.update_indices()
    return retval

class ReactionRule(object):

    def __init__(self, reactants, products, opts=[]):
        self.__reactants = reactants
        self.__products = products
        self.__options = opts

        self.initialize()

    def reactants(self):
        return copy.deepcopy(self.__reactants)

    def products(self):
        return copy.deepcopy(self.__products)

    def options(self):
        # return copy.deepcopy(self.__options)
        return self.__options

    def num_reactants(self):
        return len(self.__reactants)

    def num_products(self):
        return len(self.__products)

    def initialize(self):
        self.__correspondences = []
        self.__removed = []

        reactant_subunits = list(itertools.chain(
            *[sp.subunits for sp in self.__reactants]))
        product_subunits = list(itertools.chain(
            *[sp.subunits for sp in self.__products]))
        num_subunits = len(reactant_subunits)

        for i, su1 in enumerate(product_subunits):
            for j, su2 in enumerate(reactant_subunits):
                if (su1.name == su2.name
                    and set(su1.exclusions) == set(su2.exclusions)
                    and set(su1.modifications.keys())
                        == set(su2.modifications.keys())
                    and set(su1.domain_classes.keys())
                        == set(su2.domain_classes.keys())
                    and set(su1.domain_classes.values())
                        == set(su2.domain_classes.values())):
                    if len(self.__correspondences) > i:
                        # raise RuntimeError, "multiple correspondence found [%s]" % su1
                        print "WARN: multiple correspondence found [%s]" % su1
                    elif j in self.__correspondences:
                        print "WARN: multiple correspondence skipped [%s]" % su1
                    else:
                        self.__correspondences.append(j)
            if len(self.__correspondences) == i:
                self.__correspondences.append(num_subunits)
                num_subunits += 1

        for i in range(len(reactant_subunits)):
            if not i in self.__correspondences:
                self.__removed.append(i)

    def __generate(self, context, reactants):
        def serno(idx):
            value = context.get(label_subunit(idx))
            if value is None:
                raise RuntimeError, (
                    "no corresponding subunit found [%s]" % label_subunit(idx))

            i, stride1, stride2 = 0, 0, 0
            while i < len(self.__reactants):
                stride1 += len(self.__reactants[i].subunits)
                if idx < stride1:
                    return value.index + stride2
                stride2 += len(reactants[i].subunits)
                i += 1

            raise RuntimeError, (
                "an invalid subunit given [%s]" % label_subunit(idx))

        reactant_subunits = list(itertools.chain(
            *[sp.subunits for sp in self.__reactants]))
        product_subunits = list(itertools.chain(
            *[sp.subunits for sp in self.__products]))
        num_subunits = len(reactant_subunits)

        indices = [serno(i) for i in range(num_subunits)]

        if len(indices) != len(set(indices)):
            return None, None

        retval = concatenate_species(*reactants)
        cache_binding = {}
        new_correspondence = []

        for i, subunit in enumerate(product_subunits):
            correspondence = self.__correspondences[i]

            if correspondence >= len(reactant_subunits):
                retval.add_subunit(copy.deepcopy(subunit))
                target = retval.subunits[-1]
            else:
                target = retval.subunits[serno(correspondence)]
            new_correspondence.append(target.index)

            # for mod in subunit.exclusions:
            #     pass

            for mod, (state, binding) in subunit.modifications.items():
                if correspondence < len(reactant_subunits):
                    if mod == "":
                        raise RuntimeError, "an empty name for modification given."
                    elif mod[0] != "_":
                        label = label_domain(label_subunit(correspondence), mod)
                    elif mod in context.keys():
                        label = mod
                    else:
                        raise RuntimeError, "invalid name [%s] given." % mod

                    mod = context.get(label)
                    if mod is None:
                        raise RuntimeError, (
                            "no corresponding context found [%s]" % mod)

                value = target.modifications.get(mod)
                if value is None:
                    newstate, newbinding = "", ""
                else:
                    newstate, newbinding = value

                if state == "":
                    pass
                elif state[0] == "_":
                    if len(state) == 1:
                        pass
                    else:
                        newstate = context[state]
                else:
                    newstate = state

                if binding == "":
                    newbinding = ""
                elif binding[0] == "_":
                    if len(binding) == 1:
                        pass
                    else:
                        newbinding = context[binding] #XXX this seems problematic
                else:
                    stride = 0
                    for j, product in enumerate(self.__products):
                        stride += len(product.subunits)
                        if stride > i:
                            label = int(binding) * len(self.__products) + j
                            break

                    if label in cache_binding.keys():
                        newbinding = cache_binding[label]
                    else:
                        newbinding = str(retval.get_binding_stride() + 1)
                        cache_binding[label] = newbinding

                target.add_modification(mod, newstate, newbinding)

        removed = [serno(i) for i in self.__removed]
        removed.sort()
        for i in reversed(removed):
            su = retval.subunits.pop(i)
        # retval.update_indices()

        markers, stride = [], 0
        for sp in self.__products:
            markers.append(new_correspondence[stride])
            stride += len(sp.subunits)

        return check_connectivity(retval, markers)

    def check_options(self, reactants, products, context, corresp=None):
        for opt in self.__options:
            if isinstance(opt, options.Option):
                if not opt.check(reactants, products, context, corresp):
                    return opt.get(reactants, products, context, corresp)
            else:
                return opt
        return 0.0 # default

    def match(self, *reactants):
        contexts = None
        if len(self.__reactants) != len(reactants):
            return contexts

        for sp1, sp2 in zip(self.__reactants, reactants):
            contexts = sp1.match(sp2, contexts)
        return contexts

    def match_partial(self, idx, sp, contexts=None):
        if idx < 0 or idx >= len(self.__reactants):
            raise RuntimeError, "invalid index [%d] given." % (idx)

        if idx == 0:
            stride = 0
        else:
            stride = sum([
                self.__reactants[i].num_subunits() for i in range(idx)])

        return self.__reactants[idx].match(sp, contexts, stride)

    def generate(self, reactants, contexts=None):
        if type(reactants) not in (list, tuple):
            # just for the safety. remove this later
            raise RuntimeError, "invalid argument [%s] given." % (
                str(reactants))

        if contexts is None:
            contexts = self.match(*reactants)

        if contexts is None or len(contexts) == 0:
            return []
        # elif len(self.__products) == 0:
        #     return [()]

        retval = []
        for context in contexts:
            products, corresp = self.__generate(context, reactants)
            if products is None:
                continue

            if len(self.__options) > 0:
                opt = self.check_options(reactants, products, context, corresp)
                if opt is not None:
                    reaction = (copy.deepcopy(reactants), products, opt)
                    retval.append(reaction)
            else:
                reaction = (copy.deepcopy(reactants), products, None)
                retval.append(reaction)
        return retval

    def __str__(self):
        return "%s>%s" % (
            "+".join([str(sp) for sp in self.__reactants]),
            "+".join([str(sp) for sp in self.__products]))

    def is_degradation(self):
        return len(self.__products) == 0

    def is_synthesis(self):
        return len(self.__reactants) == 0

class Condition(object):

    def __init__(self):
        self.priority = 0

    def match(self, sp, contexts):
        return None

class DomainClassCondition(Condition):

    def __init__(self, key, name, values):
        Condition.__init__(self)

        self.key_subunit = key
        self.name = name

        self.size_values = len(values)
        self.named, self.unnamed = [], []
        for value in values:
            if value[0] == "_":
                if len(value) > 1:
                    if value in self.unnamed:
                        raise RuntimeError, "key [%s] is not unique." % value
                    self.unnamed.append(value)
            else:
                self.named.append(value)

    def predicator(self, subunit, name):
        values = subunit.domain_classes.get(self.name)
        return (values is not None) and (name in values)

    def generator(self, subunit):
        values = subunit.domain_classes.get(self.name)
        if values is None or len(values) < self.size_values:
            return None

        values = set(values)
        for value in self.named:
            if value not in values:
                return None
            else:
                values.remove(value)

        return list(itertools.permutations(values, self.size_unnamed))

    def match(self, sp, contexts):
        if (self.name == "" or
            (self.name[0] == "_" and not contexts.has_key(self.name))):
            raise RuntimeError, "[%s] not defined." % self.name

        retval = contexts

        unnamed = copy.copy(self.unnamed)
        for key in self.unnamed:
            if contexts.has_key(key):
                # raise RuntimeError, "key [%s] already assigned" % key
                unnamed.remove(key)
                retval = retval.filter2(self.predicator, self.key_subunit, key)
        self.size_unnamed = len(unnamed) # this must be immutable in generator

        retval = retval.product_any(
            self.generator, self.key_subunit, unnamed)
        return retval

class CommutativeCondition(Condition):

    def __init__(self, key, doms):
        Condition.__init__(self)

        self.key_subunit = key
        self.doms = doms

        for dom in doms:
            if dom[0] == "_":
                raise RuntimeError, (
                    "a commutative definition for [%s] not allowed" % (dom))

    def predicator(self, subunit):
        # for dom in self.doms:
        #     if dom not in subunit.commutatives.keys():
        #         return False
        return subunit.commutatives.is_commutative(*self.doms)

    def match(self, sp, contexts):
        return contexts.filter1(self.predicator, self.key_subunit)

class ModificationNameCondition(Condition):

    def __init__(self, key, mod):
        Condition.__init__(self)

        self.key_subunit = key
        self.mod = mod

    def generator(self, subunit):
        value = subunit.modifications.get(self.mod)
        if value is None:
            value = subunit.domain_classes.get(self.mod)
            if value is None:
                return []
            else:
                return list(value)
        elif self.mod in subunit.commutatives.keys():
            return subunit.commutatives.get_commutatives(self.mod)
        else:
            return [self.mod]

    def match(self, sp, contexts):
        if self.mod == "":
            raise RuntimeError, "an empty name for modification given."
        elif self.mod[0] == "_" and not contexts.has_key(self.mod):
            raise RuntimeError, "[%s] not defined." % self.mod

        if self.mod[0] != "_":
            label_key = label_domain(self.key_subunit, self.mod)

            retval = contexts.product2(
                self.generator, self.key_subunit, label_key)
            if retval is None or len(retval) == 0:
                return retval

            keys = [key for key in retval.keys()
                if is_label_domain(key, self.key_subunit) and key != label_key]
            if len(keys) > 0:
                return retval.filter_unique(label_key, keys) #XXX: an irregular use of Contexts
            else:
                return retval
        elif contexts.has_key(self.mod):
            return contexts
        else:
            raise RuntimeError, "invalid name [%s] given." % self.mod

class SubunitContainingCondition(Condition):

    def __init__(self, key, name):
        Condition.__init__(self)
        self.key_subunit = key
        self.name = name

    def predicator(self, subunit, name):
        return (subunit.name == name)

    def modifier(self, subunit):
        return subunit.name

    def match(self, sp, contexts):
        if self.name[0] == "_":
            pttrn = copy.copy(sp.subunits)
            retval = contexts.product(self.key_subunit, pttrn)
            if len(self.name) == 1 or len(pttrn) == 0:
                return retval
            elif contexts.has_key(self.name):
                return retval.filter2(self.predicator, self.key_subunit, self.name)
            else:
                return retval.update(self.modifier, self.key_subunit, self.name)
        else:
            pttrn = [subunit for subunit in sp.subunits
                if subunit.name == self.name]
            return contexts.product(self.key_subunit, pttrn)

class ModificationStateCondition(Condition):

    def __init__(self, key, mod, state):
        Condition.__init__(self)
        self.key_subunit = key
        self.mod = mod
        self.state = state

    def predicator1(self, subunit, mod):
        value = subunit.modifications.get(mod)
        return (value is not None and value[0] == self.state)

    def predicator2(self, subunit, mod):
        value = subunit.modifications.get(mod)
        return value is not None

    def predicator3(self, subunit, mod, state):
        value = subunit.modifications.get(mod)
        return value[0] == state

    def modifier(self, subunit, mod):
        return subunit.modifications.get(mod)[0]

    def match(self, sp, contexts):
        if self.mod == "":
            raise RuntimeError, "an empty name for modification given."
        elif self.mod[0] != "_":
            label = label_domain(self.key_subunit, self.mod)
        elif contexts.has_key(self.mod):
            label = self.mod
        else:
            raise RuntimeError, "invalid name [%s] given." % self.mod

        if self.state[0] != "_":
            return contexts.filter2(
                self.predicator1, self.key_subunit, label)
        elif len(self.state) == 1: # self.state == "_"
            return contexts.filter2(
                self.predicator2, self.key_subunit, label)
        elif contexts.has_key(self.state):
            return contexts.filter3(
                self.predicator3, self.key_subunit, label, self.state)
        else:
            return contexts.update2(
                self.modifier, self.key_subunit, label, self.state)

class ModificationBindingCondition(Condition):

    def __init__(self, key, mod, binding):
        Condition.__init__(self)
        self.key_subunit = key
        self.mod = mod
        self.binding = binding

        if self.binding == "":
            self.key_binding = ""
        elif self.binding[0] == "_":
            self.key_binding = self.binding
        else:
            self.key_binding = label_binding(self.binding)

    def predicator1(self, subunit, mod):
        check_binding = (
            (lambda x: x != "") if self.binding == "_" else
            (lambda x: x == ""))
        value = subunit.modifications.get(mod)
        return (value is not None and check_binding(value[1]))

    def predicator2(self, subunit, mod, target):
        value = subunit.modifications.get(mod)
        return (value is not None and value[1] == target)

    def modifier1(self, subunit, mod):
        value = subunit.modifications.get(mod)
        if (value is not None and value[1] != ""):
            return value[1]
        else:
            return None

    def modifier2(self, subunit, mod):
        value = subunit.modifications.get(mod)
        if (value is not None and value[1] != ""):
            return value[1]
        else:
            return None

    def match(self, sp, contexts):
        if self.mod == "":
            raise RuntimeError, "an empty name for modification given."
        elif self.mod[0] != "_":
            label = label_domain(self.key_subunit, self.mod)
        elif contexts.has_key(self.mod):
            label = self.mod
        else:
            raise RuntimeError, "invalid name [%s] given." % self.mod

        if self.binding == "_" or self.binding == "":
            return contexts.filter2(
                self.predicator1, self.key_subunit, label)
        elif contexts.has_key(self.key_binding):
            return contexts.filter3(
                self.predicator2, self.key_subunit, label, self.key_binding)
        else:
            return contexts.update2(
                self.modifier1, self.key_subunit, label, self.key_binding)

class ExcludedModificationCondition(Condition):

    def __init__(self, key, mod):
        Condition.__init__(self)
        self.key_subunit = key
        self.mod = mod

    def predicator(self, subunit):
        if subunit.modifications.get(self.mod) is not None:
            return False
        else:
            return (subunit.domain_classes.get(self.mod) is None)

    def match(self, sp, contexts):
        return contexts.filter1(self.predicator, self.key_subunit)

class Contexts(object):

    def __init__(self):
        self.__data = []
        self.__keys = None

    def __iter__(self):
        return iter(self.__data)

    def __len__(self):
        return len(self.__data)

    def __str__(self):
        return str(self.__data)

    def initialize(self):
        if self.__keys is not None:
            raise RuntimeError, "initialized called twice."

        self.__keys = []
        self.__data = [{}]

    def clear_locals(self):
        if self.__keys is None:
            return

        for i, key in reversed(zip(range(len(self.__keys)), self.__keys)):
            if key[0] != "_" and (len(key) > 7 and key[: 7] != "subunit"):
                for data in self.__data:
                    del data[key]
                self.__keys.pop(i)

    def _append(self, value):
        if not isinstance(value, dict):
            raise RuntimeError

        if self.__keys is None:
            self.__keys = value.keys()
        else:
            if len(value) != len(self.__keys):
                raise RuntimeError, "invalid keys [%s] given." % (value.keys())
            for key in value.keys():
                if not key in self.__keys:
                    raise RuntimeError, "invalid key [%s] found." % (key)

        self.__data.append(value)

    def has_key(self, key):
        return (key in self.__keys)

    def keys(self):
        return copy.copy(self.__keys)

    def __str__(self):
        return str(self.__data)

    def num_subunits(self):
        indices = [
            int(key[7: ]) for key in self.__keys if key[: 7] == "subunit"]
        if len(indices) == 0:
            return 0
        else:
            return max(indices) + 1
        # return len([key for key in self.__keys if key[: 7] == "subunit"])

    def product(self, key, values):
        """key is always a subunit."""
        if self.has_key(key):
            raise RuntimeError, "key [%s] already exists." % (key)

        retval = Contexts()
        for context in self.__data:
            for value in values:
                newcontext = copy.copy(context)
                newcontext[key] = value
                retval._append(newcontext)
        return retval

    def product2(self, generator, key1, key2):
        """key1 is always a subunit."""
        if self.has_key(key2):
            raise RuntimeError, "key [%s] already exists." % (key2)
        elif not self.has_key(key1):
            raise RuntimeError, "invalid key [%s] found." % (key1)

        retval = Contexts()
        for context in self.__data:
            values = generator(context[key1])
            for value in values:
                if value is None:
                    continue
                newcontext = copy.copy(context)
                newcontext[key2] = value
                retval._append(newcontext)
        return retval

    def product_any(self, generator, key, newkeys):
        if not self.has_key(key):
            raise RuntimeError, "invalid key [%s] given." % (key)
        elif any([self.has_key(elem) for elem in newkeys]):
            raise RuntimeError, "key [%s] already exists." % str(newkeys)

        retval = Contexts()
        for context in self.__data:
            subsets = generator(context[key])
            if subsets is None:
                continue
            for values in subsets:
                if values is None or len(values) != len(newkeys):
                    raise RuntimeError, "invalid return value [%s]" % str(values)
                newcontext = copy.copy(context)
                for newkey, value in zip(newkeys, values):
                    newcontext[newkey] = value
                retval._append(newcontext)
        return retval

    def filter1(self, predicator, key1):
        """key1 is always a subunit."""
        if not self.has_key(key1):
            raise RuntimeError, "invalid key [%s] found." % (key1)

        retval = Contexts()
        for context in self.__data:
            if predicator(context[key1]):
                retval._append(context)
        return retval

    def filter2(self, predicator, key1, key2):
        """key1 always indicates a subunit, but key2 doesn't."""
        if not self.has_key(key1):
            raise RuntimeError, "invalid key [%s] found." % (key1)
        if not self.has_key(key2):
            raise RuntimeError, "invalid key [%s] found." % (key2)

        retval = Contexts()
        for context in self.__data:
            if predicator(context[key1], context[key2]):
                retval._append(context)
        return retval

    def filter3(self, predicator, key1, key2, key3):
        """key1 always indicates a subunit, but key2 doesn't."""
        if not self.has_key(key1):
            raise RuntimeError, "invalid key [%s] found." % (key1)
        if not self.has_key(key2):
            raise RuntimeError, "invalid key [%s] found." % (key2)
        if not self.has_key(key3):
            raise RuntimeError, "invalid key [%s] found." % (key3)

        retval = Contexts()
        for context in self.__data:
            if predicator(context[key1], context[key2], context[key3]):
                retval._append(context)
        return retval

    def filter_unique(self, newkey, keys):
        #XXX: an irregular use of Contexts

        if not self.has_key(newkey):
            raise RuntimeError, "invalid key [%s] found." % (newkey)
        for key in keys:
            if not self.has_key(key):
                raise RuntimeError, "invalid key [%s] found." % (key)

        retval = Contexts()
        for context in self.__data:
            newvalue = context[newkey]
            if all([context[key] != newvalue for key in keys]):
                retval._append(context)
        return retval

    def update(self, modifier, key1, key2):
        """key1 always indicates a subunit, but key2 doesn't."""
        return self.product2(lambda su: [modifier(su)], key1, key2)

    def update2(self, modifier, key1, key2, key3):
        if not self.has_key(key1):
            raise RuntimeError, "invalid key [%s] found." % (key1)
        if not self.has_key(key2):
            raise RuntimeError, "invalid key [%s] found." % (key2)
        if self.has_key(key3):
            raise RuntimeError, "key [%s] already exists." % (key3)

        retval = Contexts()
        for context in self.__data:
            value = modifier(context[key1], context[key2])
            if value is not None:
                newcontext = copy.copy(context)
                newcontext[key3] = value
                retval._append(newcontext)
        return retval

class CmpSubunit:

    def __init__(self, sp):
        self.__species = sp

        self.initialize()

    def initialize(self):
        self.__species.update_indices()
        self.__subunits = copy.deepcopy(self.__species.subunits)

        self.__bindings = {}
        for i, su in enumerate(self.__subunits):
            for mod, (state, binding) in su.modifications.items():
                if binding == "":
                    continue
                elif binding[0] == "_":
                    continue

                if not binding in self.__bindings.keys():
                    self.__bindings[binding] = [(i, mod)]
                elif len(self.__bindings[binding]) == 1:
                    self.__bindings[binding].append((i, mod))
                else:
                    raise RuntimeError, "an invalid bindig found. [%s]" % (binding)

    def sort_recurse(self, idx, stride):
        su = self.__species.subunits[idx]
        if su.index < 0:
            su.index = stride
            stride += 1
        else:
            return stride

        mods = su.modifications.keys()
        for mod in sorted(mods):
            state, binding = su.modifications[mod]
            if binding != "" and binding[0] != "_":
                pair = self.__bindings[binding]
                tgt_idx, tgt_mod = pair[0] if pair[1][0] == idx else pair[1]
                stride = self.sort_recurse(tgt_idx, stride)
        return stride

    def sort(self):
        """only available for a connected graph"""
        self.__species.subunits.sort(cmp=self)
        self.initialize()

        for su in self.__species.subunits:
            su.index = -1

        self.sort_recurse(0, 0)
        self.__species.subunits.sort(key=lambda su: su.index)
        self.__species.update_indices()

        stride, newbindings = 1, {}
        for su in self.__species.subunits:
            mods = su.modifications.keys()
            for mod in sorted(mods):
                state, binding = su.modifications[mod]
                if binding == "" or binding[0] == "_":
                    continue

                newbinding = newbindings.get(binding)

                #XXX: updating subunits through the reference)
                if newbinding is None:
                    su.modifications[mod] = (state, str(stride))
                    newbindings[binding] = str(stride)
                    stride += 1
                else:
                    su.modifications[mod] = (state, newbinding)

    def cmp_recurse(self, idx1, idx2, ignore):
        if idx1 == idx2:
            return 0
        elif idx1 > idx2:
            pair_key = (idx1, idx2)
        else:
            pair_key = (idx2, idx1)

        if pair_key in ignore:
            return 0 # already checked

        su1, su2 = self.__subunits[idx1], self.__subunits[idx2]
        if su1.name != su2.name:
            return cmp(su1.name, su2.name)

        mods1, mods2 = su1.exclusions, su2.exclusions
        if len(mods1) != len(mods2):
            return cmp(len(mods1), len(mods2))

        mods1.sort()
        mods2.sort()

        for mod1, mod2 in zip(mods1, mods2):
            if mod1 != mod2:
                return cmp(mod1, mod2)

        mods1, mods2 = su1.modifications.keys(), su2.modifications.keys()
        if len(mods1) != len(mods2):
            return cmp(len(mods1), len(mods2))

        mods1.sort()
        mods2.sort()

        for mod1, mod2 in zip(mods1, mods2):
            if mod1 != mod2:
                return cmp(mod1, mod2)

            state1, state2 = (
                su1.modifications[mod1][0], su2.modifications[mod2][0])
            if state1 != state2:
                return cmp(state1, state2)

        ignore.append(pair_key)

        for mod1, mod2 in zip(mods1, mods2):
            binding1, binding2 = (
                su1.modifications[mod1][1], su2.modifications[mod2][1])
            if binding1 == binding2:
                continue
            elif binding1 == "" or binding2 == "":
                ignore.pop()
                return cmp(binding1, binding2)

            pair1, pair2 = self.__bindings[binding1], self.__bindings[binding2]
            tgt_idx1, tgt_mod1 = pair1[0] if pair1[1][0] == idx1 else pair1[1]
            tgt_idx2, tgt_mod2 = pair2[0] if pair2[1][0] == idx2 else pair2[1]

            if tgt_mod1 != tgt_mod2:
                ignore.pop()
                return cmp(tgt_mod1, tgt_mod2)

            retval = self.cmp_recurse(tgt_idx1, tgt_idx2, ignore)
            if retval != 0:
                ignore.pop()
                return retval
        else:
            ignore.pop()
            return 0

    def __call__(self, lhs, rhs):
        return self.cmp_recurse(lhs.index, rhs.index, [])


if __name__ == "__main__":
    s1 = Species()
    su1 = Subunit("X")
    su1.add_modification("a", binding=1)
    su1.add_modification("b", state="c")
    s1.add_subunit(su1)
    su2 = Subunit("X")
    su2.add_modification("a", binding=2)
    su2.add_modification("b", state="c", binding=1)
    s1.add_subunit(su2)
    su3 = Subunit("Y")
    su3.add_modification("d", binding=2)
    su3.add_modification("b", state="f")
    su3.add_modification("g")
    s1.add_subunit(su3)

    s2 = Species()
    su1 = Subunit("X")
    su1.add_modification("a", binding="_")
    su1.add_modification("b", state="c")
    s2.add_subunit(su1)

    s3 = Species()
    su1 = Subunit("X")
    su1.add_modification("a", binding=1)
    su1.add_modification("b", binding="_")
    s3.add_subunit(su1)
    su2 = Subunit("Y")
    su2.add_modification("d", binding=1)
    s3.add_subunit(su2)

    s4 = Species()
    su1 = Subunit("X")
    su1.add_modification("a", binding="_")
    # su1.add_modification("a", binding="_1")
    s4.add_subunit(su1)

    s5 = Species()
    su1 = Subunit("X")
    su1.add_modification("a")
    s5.add_subunit(su1)

    s6 = Species()
    su1 = Subunit("X")
    su1.add_modification("b", binding=1)
    s6.add_subunit(su1)
    su2 = Subunit("X")
    su2.add_modification("a", binding=1)
    s6.add_subunit(su2)

    s7 = Species()
    su1 = Subunit("_")
    su1.add_modification("b", binding="")
    s7.add_subunit(su1)

    s8 = Species()
    su1 = Subunit("_1")
    su1.add_modification("a", binding="1")
    s8.add_subunit(su1)
    su2 = Subunit("_2")
    su2.add_modification("d", binding="1")
    s8.add_subunit(su2)

    s9 = Species()
    su1 = Subunit("X")
    su1.add_modification("a", binding="_1")
    s9.add_subunit(su1)
    su2 = Subunit("Y")
    su2.add_modification("d", binding="_1")
    s9.add_subunit(su2)

    s10 = Species()
    su1 = Subunit("X")
    su1.add_modification("b", state="_1")
    s10.add_subunit(su1)

    s11 = Species()
    su1 = Subunit("X")
    su1.add_modification("b", state="_")
    s11.add_subunit(su1)

    print "Species1 =", s1

    for pttrn in [s2, s3, s4, s5, s6, s7, s8, s9, s10, s11]:
        print "Species Pattern =", pttrn
        contexts = pttrn.match(s1)
        num_patterns = len(contexts)
        print "Match Result => %d Patterns found: %s" % (num_patterns, contexts)

    # Grb2(SH2!1,SH3!2).Grb2(SH2!3,SH3!4).Grb2(SH2!5,SH3!6).Grb2(SH2!7,SH3!8).Shc(PTB!9,Y317~pY!3).Shc(PTB!10,Y317~pY!7).Sos(dom!2).Sos(dom!4).Sos(dom!6).Sos(dom!8).egf(r!11).egf(r!12).egfr(Y1068~pY!1,Y1148~pY!9,l!11,r!13).egfr(Y1068~pY!5,Y1148~pY!10,l!12,r!13)
    # egfr(Y1148~pY!1).Shc(PTB!1,Y317~pY)

    sp1 = Species()
    su = Subunit("Grb2")
    su.add_modification("SH2", binding=1)
    su.add_modification("SH3", binding=2)
    sp1.add_subunit(su)
    su = Subunit("Grb2")
    su.add_modification("SH2", binding=3)
    su.add_modification("SH3", binding=4)
    sp1.add_subunit(su)
    su = Subunit("Grb2")
    su.add_modification("SH2", binding=5)
    su.add_modification("SH3", binding=6)
    sp1.add_subunit(su)
    # su = Subunit("Grb2")
    # su.add_modification("SH2", binding=7)
    # su.add_modification("SH3", binding=8)
    # sp1.add_subunit(su)
    su = Subunit("Shc")
    su.add_modification("PTB", binding=9)
    su.add_modification("Y317", state="pY", binding=3)
    sp1.add_subunit(su)
    su = Subunit("Shc")
    su.add_modification("PTB", binding=10)
    # su.add_modification("Y317", state="pY", binding=7)
    su.add_modification("Y317", state="Y")
    sp1.add_subunit(su)
    su = Subunit("Sos")
    su.add_modification("dom", binding=2)
    sp1.add_subunit(su)
    su = Subunit("Sos")
    su.add_modification("dom", binding=4)
    sp1.add_subunit(su)
    su = Subunit("Sos")
    su.add_modification("dom", binding=6)
    sp1.add_subunit(su)
    # su = Subunit("Sos")
    # su.add_modification("dom", binding=8)
    # sp1.add_subunit(su)
    su = Subunit("egf")
    su.add_modification("r", binding=11)
    sp1.add_subunit(su)
    su = Subunit("egf")
    su.add_modification("r", binding=12)
    sp1.add_subunit(su)
    su = Subunit("egfr")
    su.add_modification("Y1068", state="pY", binding=1)
    su.add_modification("Y1148", state="pY", binding=9)
    su.add_modification("l", binding=11)
    su.add_modification("r", binding=13)
    sp1.add_subunit(su)
    su = Subunit("egfr")
    su.add_modification("Y1068", state="pY", binding=5)
    su.add_modification("Y1148", state="pY", binding=10)
    su.add_modification("l", binding=12)
    su.add_modification("r", binding=13)
    sp1.add_subunit(su)

    sp2 = Species()
    su = Subunit("egfr")
    su.add_modification("Y1068", state="pY", binding=1)
    sp2.add_subunit(su)
    su = Subunit("Grb2")
    su.add_modification("SH2", binding=1)
    su.add_modification("SH3", binding="_")
    sp2.add_subunit(su)

    print ""
    print sp1
    print sp2
    print sp2.match(sp1)

    sp3 = Species()
    su = Subunit("egfr")
    su.add_modification("r", binding="_")
    su.add_modification("Y1148", state="pY", binding=1)
    sp3.add_subunit(su)
    su = Subunit("Shc")
    su.add_modification("PTB", binding=1)
    su.add_modification("Y317", state="Y")
    sp3.add_subunit(su)

    print sp3
    print sp3.match(sp1)

    sp4 = Species()
    su = Subunit("Shc")
    su.add_modification("PTB", binding=1)
    su.add_modification("Y317", state="pY")
    sp4.add_subunit(su)
    su = Subunit("egfr")
    su.add_modification("r", binding="_")
    su.add_modification("Y1148", state="pY", binding=1)
    sp4.add_subunit(su)

    rr1 = FirstOrderReactionRule(sp3, sp4)
    print rr1.match(sp1)
