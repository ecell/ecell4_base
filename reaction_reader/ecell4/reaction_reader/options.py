class Option(object):

    def __init__(self):
        pass

    def check(self, reactants, products, context, corresp=None):
        return None

class IncludeReactants(Option):

    def __init__(self, idx, pttrn):
        Option.__init__(self)

        if type(pttrn) != str:
            raise RuntimeError
        elif pttrn[0] == "_":
            raise RuntimeError

        self.__idx = idx
        self.__pttrn = pttrn

    def check(self, reactants, products, context, corresp=None):
        if not (len(reactants) > self.__idx - 1):
            print reactants
            raise RuntimeError

        sp = reactants[self.__idx - 1]
        return (self.__pttrn in [su.name for su in sp.subunits])

class ExcludeReactants(Option):

    def __init__(self, idx, pttrn):
        Option.__init__(self)

        if type(pttrn) != str:
            raise RuntimeError
        elif pttrn[0] == "_":
            raise RuntimeError

        self.__idx = idx
        self.__pttrn = pttrn

    def check(self, reactants, products, context, corresp=None):
        if not (len(reactants) > self.__idx - 1):
            print reactants
            raise RuntimeError

        sp = reactants[self.__idx - 1]
        return not (self.__pttrn in [su.name for su in sp.subunits])

class IncludeProducts(Option):

    def __init__(self, idx, pttrn):
        Option.__init__(self)

        if type(pttrn) != str:
            raise RuntimeError
        elif pttrn[0] == "_":
            raise RuntimeError

        self.__idx = idx
        self.__pttrn = pttrn

    def check(self, reactants, products, context, corresp=None):
        if corresp is None:
            if len(products) < self.__idx:
                raise RuntimeError
            sp = products[self.__idx]
        else:
            if len(corresp) < self.__idx:
                raise RuntimeError
            sp = products[corresp[self.__idx - 1]]
        return (self.__pttrn in [su.name for su in sp.subunits])

class ExcludeProducts(Option):

    def __init__(self, idx, pttrn):
        Option.__init__(self)

        if type(pttrn) != str:
            raise RuntimeError
        elif pttrn[0] == "_":
            raise RuntimeError

        self.__idx = idx
        self.__pttrn = pttrn

    def check(self, reactants, products, context, corresp=None):
        if corresp is None:
            if len(products) < self.__idx:
                raise RuntimeError
            sp = products[self.__idx]
        else:
            if len(corresp) < self.__idx:
                raise RuntimeError
            sp = products[corresp[self.__idx - 1]]
        return not (self.__pttrn in [su.name for su in sp.subunits])
