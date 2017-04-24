from __future__ import print_function
from collections import defaultdict

import numbers

from ecell4 import Species, ReactionRule
from ecell4.ode import ODEReactionRule

from ecell4.util import parseobj
from ecell4.util.decorator_base import just_parse

def escape_serial(val):
    if isinstance(val, Species):
        return escape_serial(val.serial())
    return "\\left[{}\\right]".format(val)

def function_parseobj(obj):
    reserved_vars = ['_t', 'pi']
    reserved_funcs = ['exp', 'log', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan'] # , 'pow']

    if not isinstance(obj, parseobj.ParseObj) or obj._size() != 1:
        return None

    last = obj._last()
    if last.name not in reserved_funcs or last.args is None or len(last.args) != 1:
        return None

    return "{{\\{}}} {}".format(last.name, wrap_parseobj(last.args[0], True))

def wrap_parseobj(obj, force=False):
    res = function_parseobj(obj)
    if res is not None:
        return res

    if (isinstance(obj, (parseobj.AnyCallable, parseobj.ParseObj))
            or (not force and isinstance(obj, (parseobj.DivExp, parseobj.MulExp, parseobj.PowExp)))
            or (isinstance(obj, numbers.Number) and obj >= 0)
            or not isinstance(obj, (numbers.Number, parseobj.ExpBase))):
        return convert_parseobj(obj)
    return "\\left({}\\right)".format(convert_parseobj(obj))

def convert_parseobj(obj):
    res = function_parseobj(obj)
    if res is not None:
        return res

    if isinstance(obj, (parseobj.AnyCallable, parseobj.ParseObj)):
        return escape_serial(str(obj))
    elif isinstance(obj, parseobj.PosExp):
        return "+{}".format(wrap_parseobj(obj._target))
    elif isinstance(obj, parseobj.NegExp):
        return "-{}".format(wrap_parseobj(obj._target))
    elif isinstance(obj, parseobj.AddExp):
        return "+".join([convert_parseobj(elem) for elem in obj._elements()])
    elif isinstance(obj, parseobj.SubExp):
        elems = obj._elements()
        res = convert_parseobj(elems[0])
        res += "".join(["-{}".format(wrap_parseobj(elem)) for elem in elems[1: ]])
        return res
    elif isinstance(obj, parseobj.DivExp):
        elems = obj._elements()
        res = convert_parseobj(elems[0])
        for elem in elems[1: ]:
            res = "\\frac{{{}}}{{{}}}".format(res, convert_parseobj(elem))
        return res
    elif isinstance(obj, parseobj.MulExp):
        return "".join([wrap_parseobj(elem) for elem in obj._elements()])
    elif isinstance(obj, parseobj.PowExp):
        return "{{{}}}^{{{}}}".format(wrap_parseobj(obj._lhs, True), convert_parseobj(obj._rhs))
    elif isinstance(obj, parseobj.ExpBase):
        raise ValueError("The given type of ParseObj is not supported yet [{}].".format(repr(obj)))

    return str(obj)

def equations(m):
    derivatives = defaultdict(list)
    equations = {}
    for i, rr in enumerate(m.reaction_rules()):
        name = "v_{{{}}}".format(i + 1)

        if isinstance(rr, ReactionRule) or (isinstance(rr, ODEReactionRule) and rr.is_massaction()):
            equations[name] = "{}{}".format(rr.k(), "".join([escape_serial(sp) for sp in rr.reactants()]))
        elif isinstance(rr, ODEReactionRule):
            #XXX: rr.is_massaction() is False
            rl = rr.get_ratelaw()
            parser = just_parse()
            parser.set_callback()
            obj = parser.eval(rl.as_string(), params={'pow': pow})
            # obj = just_parse()._ParseDecorator__evaluate(rl.as_string())
            equations[name] = convert_parseobj(obj)

        stoich = defaultdict(int)

        if isinstance(rr, ReactionRule):
            for sp in rr.reactants():
                stoich[escape_serial(sp)] -= 1
            for sp in rr.products():
                stoich[escape_serial(sp)] += 1

            for serial, coef in stoich.items():
                if coef == 0:
                    continue
                elif coef == 1:
                    pref = ""
                elif coef == -1:
                    pref = "-"
                else: # coef.is_integer():
                    pref = "{:d} ".format(coef)

                if len(derivatives[serial]) == 0 and coef > 0:
                    derivatives[serial].append("{}{}".format(pref, name))
                elif coef > 0:
                    derivatives[serial].append("+{}{}".format(pref, name))
                else:
                    derivatives[serial].append("{}{}".format(pref, name))
        elif isinstance(rr, ODEReactionRule):
            for sp, coef in zip(rr.reactants(), rr.reactants_coefficients()):
                stoich[escape_serial(sp)] -= coef
            for sp, coef in zip(rr.products(), rr.products_coefficients()):
                stoich[escape_serial(sp)] += coef

            for serial, coef in stoich.items():
                if coef == 0:
                    continue
                elif coef == 1:
                    pref = ""
                elif coef == -1:
                    pref = "-"
                elif coef.is_integer():
                    pref = "{:d} ".format(int(coef))
                else:
                    pref = "{:g} ".format(coef)

                if len(derivatives[serial]) == 0 and coef > 0:
                    derivatives[serial].append("{}{}".format(pref, name))
                elif coef > 0:
                    derivatives[serial].append("+{}{}".format(pref, name))
                else:
                    derivatives[serial].append("{}{}".format(pref, name))

    res = []
    for serial in sorted(derivatives.keys()):
        res.append(
            "\\frac{{\\mathrm{{d}} {}}}{{\\mathrm{{d}} t}} &=& {}".format(serial, ''.join(derivatives[serial])))
    for name in sorted(equations.keys()):
        res.append(
            "{} &=& {}".format(name, equations[name]))

    res = "\\begin{{eqnarray}}\n{}\n\\end{{eqnarray}}".format(',\\\\\n'.join(res))
    return res


if __name__ == "__main__":
    from ecell4 import reaction_rules, get_model

    with reaction_rules():
        A + B == C | (1.0, 2.0)
        C > D | 3.0
        D > A + B | 4.0

    m = get_model()

    print(equations(m))
