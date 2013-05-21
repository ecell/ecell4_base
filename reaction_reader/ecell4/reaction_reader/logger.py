from functools import wraps


def log_call(func):
    bits = []
    for attr in ['__module__', '__name__']:
        if hasattr(func, attr):
            bits.append(getattr(func, attr))
    func_id = '.'.join(bits)
    @wraps(func)
    def wrapped(*args, **kwargs):
        print "%s, args=%s, kwargs=%s" % (func_id, args, kwargs) # XXX: use logging
        ret = func(*args, **kwargs)
        print "%s, retval=%s" % (func_id, ret) # XXX: use logging
        return ret
    return wrapped
