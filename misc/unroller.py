import re
import sys

sys.argv.pop(0)

def cb(g):
    ind = g.group(1)
    exc = g.group(2)
    cond = g.group(3)
    vars = []
    for m in re.finditer(r'(?:^|[^a-zA-Z0-9_])((?:[a-zA-Z_][a-zA-Z0-9_]*|->)+)', cond):
        var = m.group(1)
        if var not in vars:
            vars.append(var)
    vars_fmt = []
    for var in vars:
        vars_fmt.append('%s=%%g' % var)
    vars = ' % '.join(vars)
    vars_fmt = ', '.join(vars_fmt)
    return '''%(ind)sif (!(%(cond)s))
%(ind)s{
%(ind)s    throw %(exc)s((boost::format("%(cond)s: %(vars_fmt)s") %% %(vars)s).str());
%(ind)s}
''' % locals()

for f in sys.argv:
    c = re.sub(r'''([\t]*)THROW_UNLESS\(\s*([^, ]*)\s*,\s*([^,)]*)\s*\);''', cb,
            file(f).read())
    file(f, 'w').write(c)
