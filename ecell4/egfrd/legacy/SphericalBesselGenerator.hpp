#ifndef __SPHERICALBESSELGENERATOR_HPP
#define __SPHERICALBESSELGENERATOR_HPP

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>

#include "Defs.hpp"

// #define NO_BESSEL_TABLE

#ifdef NO_BESSEL_TABLE
#include "tablegen/sjy_table.hpp"

namespace sb_table
{

struct Table
{
    unsigned int N;
    double x_start;
    double delta_x;
    std::vector<double> y;
};

const unsigned int sj_table_min = 4;
const unsigned int sj_table_max = 51;
const unsigned int sy_table_min = 3;
const unsigned int sy_table_max = 40;
const unsigned int sjy_table_resolution = 35;

} // sb_table
#endif

class SphericalBesselGenerator
{

    typedef UnsignedInteger Index;

public:

    SphericalBesselGenerator()
    {
#ifdef NO_BESSEL_TABLE
        // std::cout << "SphericalBesselGenerator::SphericalBesselGenerator() was called."<< std::endl;
        sjy_table table = jnyn(std::max(sb_table::sj_table_max, sb_table::sy_table_max), sb_table::sjy_table_resolution);

        sj_table_.resize(sb_table::sj_table_max + 1);
        for (unsigned int n(sb_table::sj_table_min); n<= sb_table::sj_table_max; ++n)
        {
            const int start(searchsorted(table.z, minz_j(n)));
            const double z_start(table.z.at(start));
            const int end(searchsorted(table.z, maxz_j(n)));
            const std::vector<double> js(get_sub_sequence_from_matrix2(table.j, table.jdot, n, start, end));

            const sb_table::Table sj_table_n = {end - start, z_start, table.delta, js};
            sj_table_[n] = sj_table_n;
        }

        sy_table_.resize(sb_table::sy_table_max + 1);
        for (unsigned int n(sb_table::sy_table_min); n<= sb_table::sy_table_max; ++n)
        {
            const int start(searchsorted(table.z, minz_y(n)));
            const double z_start(table.z.at(start));
            const int end(searchsorted(table.z, maxz_y(n)));
            const std::vector<double> ys(get_sub_sequence_from_matrix2(table.y, table.ydot, n, start, end));

            const sb_table::Table sy_table_n = {end - start, z_start, table.delta, ys};
            sy_table_[n] = sy_table_n;
        }
        // std::cout << "SphericalBesselGenerator::SphericalBesselGenerator() was done."<< std::endl;
#endif
    }

    ~SphericalBesselGenerator()
    {
        ; // do nothing
    }

    Real j(UnsignedInteger n, Real z) const;

    Real y(UnsignedInteger n, Real z) const;

    static UnsignedInteger getMinNJ();
    static UnsignedInteger getMinNY();
    static UnsignedInteger getMaxNJ();
    static UnsignedInteger getMaxNY();

    static SphericalBesselGenerator const& instance();

#ifdef NO_BESSEL_TABLE
    Real _j_table(UnsignedInteger n, Real z) const;
    Real _y_table(UnsignedInteger n, Real z) const;
    sb_table::Table const* getSJTable(UnsignedInteger n) const;
    sb_table::Table const* getSYTable(UnsignedInteger n) const;

private:

    std::vector<sb_table::Table> sj_table_;
    std::vector<sb_table::Table> sy_table_;
#endif
};




#endif /* __SPHERICALBESSELGENERATOR_HPP */
