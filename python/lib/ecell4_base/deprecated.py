import warnings
import functools

def deprecated(suggest=None):
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            warnings.simplefilter('always', DeprecationWarning)
            warnings.warn("{} is deprecated.".format(func.__name__),
                          category = DeprecationWarning,
                          stacklevel = 2)
            warnings.simplefilter('default', DeprecationWarning)
            return func(*args, **kwargs)

        if suggest is None:
            doc = "[Deprecated]\n"
        else:
            doc = "[Deprecated] Use '" + suggest + "' instead.\n"

        if wrapper.__doc__ is None:
            wrapper.__doc__ = doc
        else:
            wrapper.__doc__ = doc + wrapper.__doc__

        return wrapper
    return decorator

