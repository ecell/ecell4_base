#ifndef __CYLINDRICALBESSELGENERATOR_HPP
#define __CYLINDRICALBESSELGENERATOR_HPP

#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#include "Defs.hpp"

#define NO_BESSEL_TABLE

#ifdef NO_BESSEL_TABLE
#include "tablegen/cjy_table.hpp"

namespace cb_table
{

struct Table
{
    unsigned int N;
    double x_start;
    double delta_x;
    std::vector<double> y;
};

const unsigned int cj_table_min = 0;
const unsigned int cj_table_max = 50;
const unsigned int cy_table_min = 0;
const unsigned int cy_table_max = 50;
const unsigned int cjy_table_resolution = 35;

} // sb_table
#endif


class CylindricalBesselGenerator
{

    typedef UnsignedInteger Index;

public:

    CylindricalBesselGenerator()
    {
#ifdef NO_BESSEL_TABLE
        // std::cout << "CylindricalBesselGenerator::CylindricalBesselGenerator() was called."<< std::endl;
        cjy_table table = JnYn(std::max(cb_table::cj_table_max, cb_table::cy_table_max), cb_table::cjy_table_resolution);

        cj_table_.reserve(cb_table::cj_table_max - cb_table::cj_table_min + 1);
        for (unsigned int n(cb_table::cj_table_min); n<= cb_table::cj_table_max; ++n)
        {
            const int start(searchsorted(table.z, minz_j(n)));
            const double z_start(table.z.at(start));
            const int end(searchsorted(table.z, maxz_j(n)));
            const std::vector<double> js(get_sub_sequence_from_matrix2(table.j, table.jdot, n, start, end));

            const cb_table::Table cj_table_n = {end - start, z_start, table.delta, js};
            cj_table_.push_back(cj_table_n);
        }

        cj_table_.reserve(cb_table::cy_table_max - cb_table::cy_table_min + 1);
        for (unsigned int n(cb_table::cy_table_min); n<= cb_table::cy_table_max; ++n)
        {
            const int start(searchsorted(table.z, minz_y(n)));
            const double z_start(table.z.at(start));
            const int end(searchsorted(table.z, maxz_y(n)));
            const std::vector<double> ys(get_sub_sequence_from_matrix2(table.y, table.ydot, n, start, end));

            const cb_table::Table cy_table_n = {end - start, z_start, table.delta, ys};
            cy_table_.push_back(cy_table_n);
        }
        // std::cout << "SphericalBesselGenerator::SphericalBesselGenerator() was done."<< std::endl;
#endif
    }

    ~CylindricalBesselGenerator()
    {
        ; // do nothing
    }

    Real J(UnsignedInteger n, Real z) const;

    Real Y(UnsignedInteger n, Real z) const;

    static UnsignedInteger getMinNJ();
    static UnsignedInteger getMinNY();
    static UnsignedInteger getMaxNJ();
    static UnsignedInteger getMaxNY();

    static CylindricalBesselGenerator const& instance();

#ifdef NO_BESSEL_TABLE
    Real _J_table(UnsignedInteger n, Real z) const;
    Real _Y_table(UnsignedInteger n, Real z) const;
    cb_table::Table const* getCJTable(UnsignedInteger n) const;
    cb_table::Table const* getCYTable(UnsignedInteger n) const;

private:

    std::vector<cb_table::Table> cj_table_;
    std::vector<cb_table::Table> cy_table_;
#endif
};




#endif /* __CYLINDRICALBESSELGENERATOR_HPP */
