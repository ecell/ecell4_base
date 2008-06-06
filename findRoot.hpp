#include <gsl/gsl_roots.h>

#include "Defs.hpp"


const Real 
findRoot( gsl_function& F,
          gsl_root_fsolver* solver,
          const Real low,
          const Real high,
          const Real tolerance,
          std::string funcName );


