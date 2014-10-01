#ifndef FIND_ROOT_HPP
#define FIND_ROOT_HPP

#include <gsl/gsl_roots.h>

#include "Defs.hpp"


Real findRoot(gsl_function const& F, gsl_root_fsolver* solver, Real low,
              Real high, Real tol_abs, Real tol_rel, char const* funcName);

#endif /* FIND_ROOT_HPP */
