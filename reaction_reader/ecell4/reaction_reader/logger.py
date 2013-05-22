from functools import wraps
import logging

# logging.basicConfig(level=logging.DEBUG)


def log_call(func):
    bits = []
    for attr in ['__module__', '__name__']:
        if hasattr(func, attr):
            bits.append(getattr(func, attr))
    func_id = '.'.join(bits)
    @wraps(func)
    def wrapped(*args, **kwargs):
        logging.debug("%s, args=%s, kwargs=%s" % (func_id, args, kwargs))
        ret = func(*args, **kwargs)
        logging.debug("%s, retval=%s" % (func_id, repr(ret)))
        return ret
    return wrapped
