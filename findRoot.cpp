#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdexcept>
#include <gsl/gsl_errno.h>

#include "Logger.hpp"
#include "findRoot.hpp"

static Logger& _log(Logger::get_logger("findRoot"));

Real findRoot(gsl_function const& F, gsl_root_fsolver* solver, Real low,
              Real high, Real tol_abs, Real tol_rel, char const* funcName)
{
    Real l(low);
    Real h(high);

    gsl_root_fsolver_set(solver, const_cast<gsl_function*>(&F), l, h);

    const unsigned int maxIter(100);

    unsigned int i(0);
    for (;;)
    {
        gsl_root_fsolver_iterate(solver);
        l = gsl_root_fsolver_x_lower(solver);
        h = gsl_root_fsolver_x_upper(solver);

        const int status(gsl_root_test_interval(l, h, tol_abs,
                                                  tol_rel));

        if (status == GSL_CONTINUE)
        {
            if (i >= maxIter)
            {
                gsl_root_fsolver_free(solver);
                _log.error("%s: failed to converge.", funcName);
                throw std::exception();
            }
        }
        else
        {
            break;
        }

        ++i;
    }
  

    const Real root(gsl_root_fsolver_root(solver));

    return root;
}
