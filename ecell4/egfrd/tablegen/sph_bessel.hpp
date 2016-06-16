#ifndef __ECELL4_SPH_BESSEL_HPP
#define __ECELL4_SPH_BESSEL_HPP

#include <utility>
#include <vector>
#include <limits>

typedef std::pair<std::vector<double>, std::vector<double> > values;

const double inf = std::numeric_limits<double>::infinity();

values sphj_array(const int n, const double x);
values sphy_array(const int n, const double x);

#endif
