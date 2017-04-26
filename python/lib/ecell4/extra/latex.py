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

def function_parseobj(obj, params):
    # reserved_vars = ['_t', 'pi']
    reserved_funcs = ['exp', 'log', 'sin', 'cos', 'tan', 'asin', 'acos', 'atan'] # , 'pow']

    if not isinstance(obj, parseobj.ParseObj) or obj._size() != 1:
        return None

    last = obj._last()
    if last.name not in reserved_funcs or last.args is None or len(last.args) != 1:
        return None

    return "{{\\{}}} {}".format(last.name, wrap_parseobj(last.args[0], params, True))

def wrap_parseobj(obj, params, force=False):
    res = function_parseobj(obj, params)
    if res is not None:
        return res

    if (isinstance(obj, (parseobj.AnyCallable, parseobj.ParseObj))
            or (not force and isinstance(obj, (parseobj.DivExp, parseobj.MulExp, parseobj.PowExp)))
            or (isinstance(obj, numbers.Number) and (obj >= 0 or params is not None))
            or not isinstance(obj, (numbers.Number, parseobj.ExpBase))):
        return convert_parseobj(obj, params)
    return "\\left({}\\right)".format(convert_parseobj(obj, params))

def convert_parseobj(obj, params):
    res = function_parseobj(obj, params)
    if res is not None:
        return res

    if isinstance(obj, (parseobj.AnyCallable, parseobj.ParseObj)):
        return escape_serial(str(obj))
    elif isinstance(obj, parseobj.PosExp):
        return "+{}".format(wrap_parseobj(obj._target, params))
    elif isinstance(obj, parseobj.NegExp):
        return "-{}".format(wrap_parseobj(obj._target, params))
    elif isinstance(obj, parseobj.AddExp):
        return "+".join([convert_parseobj(elem, params) for elem in obj._elements()])
    elif isinstance(obj, parseobj.SubExp):
        elems = obj._elements()
        res = convert_parseobj(elems[0], params)
        res += "".join(["-{}".format(wrap_parseobj(elem, params)) for elem in elems[1: ]])
        return res
    elif isinstance(obj, parseobj.DivExp):
        elems = obj._elements()
        res = convert_parseobj(elems[0], params)
        for elem in elems[1: ]:
            res = "\\frac{{{}}}{{{}}}".format(res, convert_parseobj(elem, params))
        return res
    elif isinstance(obj, parseobj.MulExp):
        return "".join([wrap_parseobj(elem, params) for elem in obj._elements()])
    elif isinstance(obj, parseobj.PowExp):
        return "{{{}}}^{{{}}}".format(wrap_parseobj(obj._lhs, params, True), convert_parseobj(obj._rhs, params))
    elif isinstance(obj, parseobj.ExpBase):
        raise ValueError("The given type of ParseObj is not supported yet [{}].".format(repr(obj)))

    if params is not None and isinstance(obj, numbers.Number):
        idx = len(params) + 1
        k = "k_{{{}}}".format(idx)
        params[k] = (idx, obj)
        return k
    return str(obj)

def equations(m, inline=False, constants=True):
    derivatives = defaultdict(list)
    equations = {}
    params = {} if constants else None
    for i, rr in enumerate(m.reaction_rules()):
        name = "v_{{{}}}".format(i + 1)

        if isinstance(rr, ReactionRule) or (isinstance(rr, ODEReactionRule) and rr.is_massaction()):
            if constants:
                idx = len(params) + 1
                k = "k_{{{}}}".format(idx)
                params[k] = (idx, rr.k())
            elif rr.k() != 1:
                k = "{:g}".format(rr.k())
            else:
                k = ""
            equations[name] = (i + 1, "{}{}".format(k, "".join([escape_serial(sp) for sp in rr.reactants()])))
        elif isinstance(rr, ODEReactionRule):
            #XXX: rr.is_massaction() is False
            rl = rr.get_ratelaw()
            parser = just_parse()
            parser.set_callback()
            obj = parser.eval(rl.as_string(), params={'pow': pow})
            # obj = just_parse()._ParseDecorator__evaluate(rl.as_string())
            if inline:
                equations[name] = (i + 1, wrap_parseobj(obj, params))
            else:
                equations[name] = (i + 1, convert_parseobj(obj, params))

        if inline:
            if name not in equations:
                raise ValueError("No equation is defined [{}]".format(rr.as_string()))
            _, eq = equations[name]
        else:
            eq = name

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
                    derivatives[serial].append("{}{}".format(pref, eq))
                elif coef > 0:
                    derivatives[serial].append("+{}{}".format(pref, eq))
                else:
                    # coef < 0
                    derivatives[serial].append("{}{}".format(pref, eq))
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
                    derivatives[serial].append("{}{}".format(pref, eq))
                elif coef > 0:
                    derivatives[serial].append("+{}{}".format(pref, eq))
                else:
                    # coef < 0
                    derivatives[serial].append("{}{}".format(pref, eq))

    res = []
    for serial in sorted(derivatives.keys()):
        res.append(
            "\\frac{{\\mathrm{{d}} {}}}{{\\mathrm{{d}} t}} &=& {}".format(serial, ''.join(derivatives[serial])))

    if not inline:
        for name, (_, value) in sorted(equations.items(), key=lambda x: x[1][0]):
            res.append(
                "{} &=& {}".format(name, value))

    if params:
        for name, (_, value) in sorted(params.items(), key=lambda x: x[1][0]):
            res.append(
                "{} &=& {}".format(name, value))

    res = "\\begin{{eqnarray}}\n{}\n\\end{{eqnarray}}".format(',\\\\\n'.join(res))
    return res

def display_equations(m, inline=False, constants=True):
    from IPython.display import display, Latex #, Math
    fmt = equations(m, inline, constants)
    display(Latex(fmt))


# if __name__ == "__main__":
#     from ecell4 import reaction_rules, get_model

#     with reaction_rules():
#         A + B == C | (1.0, 2.0)
#         C > D | 3.0
#         D > A + B | 4.0

#     m = get_model()

#     print(equations(m))
