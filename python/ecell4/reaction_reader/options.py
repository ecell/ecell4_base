class Option(object):

    def __init__(self):
        pass

    def check(self, reactants, products, context, corresp=None):
        return None

    def get(self, reactants, products, context, corresp=None):
        return None

class CaseIf(Option):

    def __init__(self, value, **kwargs):
        self.value = value
        self.kwargs = kwargs

        for key in self.kwargs.keys():
            if type(key) != str:
                raise RuntimeError, (
                    "a key [%s] must be a string." % (str(key)))

    def check(self, reactants, products, context, corresp=None):
        for key, value in self.kwargs.items():
            if context.get(key) != value:
                return True
        return False

    def get(self, reactants, products, context, corresp=None):
        return self.value

class CountSubunits(Option):

    def __init__(self, idx, pttrn):
        Option.__init__(self)

        if type(pttrn) != str:
            raise RuntimeError, (
                "a pattern for subunits must be a string [%s]." % (str(pttrn)))
        elif pttrn[0] == "_":
            raise RuntimeError, ("an invalid pattern [%s] given." % (pttrn))

        self.idx = idx
        self.pttrn = pttrn

class IncludeReactants(CountSubunits):

    def __init__(self, idx, pttrn):
        CountSubunits.__init__(self, idx, pttrn)

    def check(self, reactants, products, context, corresp=None):
        if not (len(reactants) > self.idx - 1):
            raise RuntimeError, (
                "the number of reactants is too small [%d < %d]." % (
                    len(reactants), self.idx))

        sp = reactants[self.idx - 1]
        return (self.pttrn in [su.name for su in sp.subunits])

class ExcludeReactants(CountSubunits):

    def __init__(self, idx, pttrn):
        CountSubunits.__init__(self, idx, pttrn)

    def check(self, reactants, products, context, corresp=None):
        if not (len(reactants) > self.idx - 1):
            raise RuntimeError, (
                "the number of reactants is too small [%d < %d]." % (
                    len(reactants), self.idx))

        sp = reactants[self.idx - 1]
        return not (self.pttrn in [su.name for su in sp.subunits])

class IncludeProducts(CountSubunits):

    def __init__(self, idx, pttrn):
        CountSubunits.__init__(self, idx, pttrn)

    def check(self, reactants, products, context, corresp=None):
        if corresp is None:
            if not (len(products) > self.idx - 1):
                raise RuntimeError, (
                    "the number of products is too small [%d < %d]." % (
                        len(products), self.idx))
            sp = products[self.idx]
        else:
            if not (len(corresp) > self.idx - 1):
                raise RuntimeError, (
                    "no corresponding subunit found [%d < %d]." % (
                        len(corresp), self.idx))
            sp = products[corresp[self.idx - 1]]
        return (self.pttrn in [su.name for su in sp.subunits])

class ExcludeProducts(CountSubunits):

    def __init__(self, idx, pttrn):
        CountSubunits.__init__(self, idx, pttrn)

    def check(self, reactants, products, context, corresp=None):
        if corresp is None:
            if not (len(products) > self.idx - 1):
                raise RuntimeError, (
                    "the number of products is too small [%d < %d]." % (
                        len(products), self.idx))
            sp = products[self.idx]
        else:
            if not (len(corresp) > self.idx - 1):
                raise RuntimeError, (
                    "no corresponding subunit found [%d < %d]." % (
                        len(corresp), self.idx))
            sp = products[corresp[self.idx - 1]]
        return not (self.pttrn in [su.name for su in sp.subunits])
