#ifndef GSL_RNG_BASE_HPP
#define GSL_RNG_BASE_HPP

#include <gsl/gsl_rng.h>

template<typename Tderived_>
struct gsl_rng_base: public gsl_rng
{
    static double _get_double(void *state)
    {
        return reinterpret_cast<Tderived_*>(state)->get_double();
    }

    static unsigned long int _get(void *state)
    {
        return reinterpret_cast<Tderived_*>(state)->get();
    }

    static void _set(void *state, unsigned long int seed)
    {
        reinterpret_cast<Tderived_*>(state)->set(seed);
    }

    gsl_rng_base()
    {
        type = &rng_type;
        state = this;
    }

    static gsl_rng_type rng_type;
};

template<typename Tderived_>
gsl_rng_type gsl_rng_base<Tderived_>::rng_type = {
    Tderived_::name,
    Tderived_::max,
    Tderived_::min,
    sizeof(Tderived_),
    &gsl_rng_base::_set,
    &gsl_rng_base::_get,
    &gsl_rng_base::_get_double
};

#endif /* GSL_RNG_BASE_HPP */
