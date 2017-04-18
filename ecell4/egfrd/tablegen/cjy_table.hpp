#ifndef EGFRD_CJY_TABLE_HPP
#define EGFRD_CJY_TABLE_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#if (_MSC_VER >= 1500)
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_nan.h>

// #include "make_table_util.hpp"
typedef std::vector<double> table_type;
typedef std::vector<std::pair<double, double> > valdot_type;

using namespace boost::numeric;

typedef ublas::matrix<double> matrix;

struct cjy_table
{
    cjy_table(table_type z, double delta,
            matrix j, matrix jdot,
            matrix y, matrix ydot) :
        z(z), delta(delta), j(j), jdot(jdot), y(y), ydot(ydot)
    {
    }

    table_type z;
    double delta;
    matrix j;
    matrix jdot;
    matrix y;
    matrix ydot;
};

inline double minz_j(const int n)
{
    // We can start table interpolation from zero because there is
    // no singularity in bessel_j for z>=0.
    return 0;
}

inline double minz_y(const int n)
{
    // return max(3., n)
    return 5.;
}

inline double maxz_j(const int n)
{
    double z = (n * n + 1) / 1.221e-4;
    if (z >= 1000)
        z = std::max(1000, n * n);
    return z;
}

inline double maxz_y(const int n)
{
    // from gsl/special/bessel_y.c:
    //  else if(GSL_ROOT3_DBL_EPSILON * x > (l*l + l + 1.0)) {
    //     int status = gsl_sf_bessel_Ynu_asympx_e(l + 0.5, x, result);
    //     ...
    double z = (n * n + n + 1) / 6.06e-6;
    // ... but this is usually too big.
    if (z >= 2000)
        z = std::max(2000, n * n);
    return z;
}

inline const table_type get_z_table(const double from, const double to, const double delta)
{
    table_type table;
    table.reserve(int((to - from) / delta));
    for (double z(from); z < to; z += delta)
        table.push_back(z);
    return table;
}

inline matrix zeros(const int length0, const int length1)
{
    matrix retval(length0, length1);
    for (int i(0); i < length0; i++)
        for (int j(0); j < length1; j++)
            retval(i, j) = 0;
    return retval;
}

inline void calculate_Jns(const int n, const double z,
        double j[], double jdot[])
{
    double *jp = new double[n+2];
    gsl_sf_bessel_Jn_array(0, n+1, z, jp);
    j[0] = jp[0];
    jdot[0] = -jp[1];
    for (int l(1); l < n+1; ++l)
    {
        j[l] = jp[l];
        jdot[l] = (jp[l-1] - jp[l+1])/2;
        //jdot[l] = (l*jp[l-1] - (l+1)*jp[l+1])/(2*l+1);
    }
    delete[] jp;
}

inline void calculate_Yns(const int n, const double z,
        double y[], double ydot[])
{
    double *yp = new double[n+2];
    gsl_sf_bessel_Yn_array(0, n+1, z, yp);
    y[0] = yp[0];
    ydot[0] = -yp[1];
    for (int l(1); l < n+1; ++l)
    {
        y[l] = yp[l];
        ydot[l] = (yp[l-1] - yp[l+1])/2;
        //ydot[l] = (n*yp[l-1] - (n+1)*yp[l])/(2*n+1);
    }
    delete[] yp;
}

inline void set_matrix(matrix& mat, const int i, const double array[], const int n)
{
    for (int index(0); index < n; ++index)
    {
        mat(i, index) = array[index];
    }
}

inline cjy_table JnYn(const int n, const int resolution)
{
    const double delta(M_PI / resolution);
    //const double delta(0.089759);
    const table_type z_table(get_z_table(
                std::min(minz_j(n), minz_y(n)),
                std::max(maxz_j(n), maxz_y(n)),
                delta));

    const int table_size(z_table.size());
    matrix j_table(zeros(table_size, n+1));
    matrix jdot_table(zeros(table_size, n+1));
    matrix y_table(zeros(table_size, n+1));
    matrix ydot_table(zeros(table_size, n+1));

    for (int i(0); i < z_table.size(); ++i)
    {
        const double z(z_table.at(i));

        // j
        double *j = new double[n+1];
        double *jdot = new double[n+1];
        calculate_Jns(n, z, j, jdot);
        set_matrix(j_table, i, j, n+1);
        set_matrix(jdot_table, i, jdot, n+1);
        delete[] j;
        delete[] jdot;

        // y
        double *y = new double[n+1];
        double *ydot = new double[n+1];
        if (z <= 0)
            for (int i(0); i < n+1; ++i)
                y[i] = ydot[i] = GSL_NEGINF;
        else
            calculate_Yns(n, z, y, ydot);
        set_matrix(y_table, i, y, n+1);
        set_matrix(ydot_table, i, ydot, n+1);
        delete[] y;
        delete[] ydot;
    }

    j_table = trans(j_table);
    jdot_table = trans(jdot_table);
    y_table = trans(y_table);
    ydot_table = trans(ydot_table);

    return cjy_table(z_table, delta, j_table, jdot_table,
            y_table, ydot_table);
}

inline int searchsorted(const table_type& z_table, const double value)
{
    int i(0);
    for (table_type::const_iterator itr(z_table.begin());
            itr != z_table.end(); ++itr)
    {
        if (value <= *itr)
            return i;
        ++i;
    }
    return i;
}

inline valdot_type get_sub_sequence_from_matrix(const matrix mat0,
        const matrix mat1, const int index0, const int start, const int end)
{
    if (start > end)
        throw "start should be no more than end.";
    valdot_type retval;
    retval.reserve(end-start);
    for (int i(start); i < end; ++i)
        retval.push_back(std::pair<double, double>(mat0(index0, i), mat1(index0, i)));
    return retval;
}

inline std::vector<double> get_sub_sequence_from_matrix2(const matrix mat0,
        const matrix mat1, const int index0, const int start, const int end)
{
    if (start > end)
        throw "start should be no more than end.";
    std::vector<double> retval;
    retval.reserve(end-start * 2);
    for (int i(start); i < end; ++i)
    {
        // retval.push_back(std::pair<double, double>(mat0(index0, i), mat1(index0, i)));
        retval.push_back(mat0(index0, i));
        retval.push_back(mat1(index0, i));
    }
    return retval;
}

#endif /* EGFRD_CJY_TABLE_HPP */
