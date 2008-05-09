#ifndef _SPHERE_HPP
#define _SPHERE_HPP

#include "sphere.hpp"

sphere<double>* __impl_sphere_new(
        const double& x, const double& y,
        const double& z, const double& r)
{
    return new sphere<double>(position<double>(x, y, z), r);
}

sphere<double>* __impl_sphere_new(
        const double (&p)[3], const double& r)
{
    return new sphere<double>(position<double>(p), r);
}


sphere<double> *__impl_sphere_clone(const sphere<double>* src)
{
    return new sphere<double>(*src);
}

#endif /* _SPHERE_HPP */
