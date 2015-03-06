#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_nan.h>
#include "make_table_util.hpp"

using namespace boost::numeric;

typedef ublas::matrix<double> matrix;


struct sjy_table
{
    sjy_table(table_type z, double delta,
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

double minz_j(const int n)
{
    // We can start table interpolation from zero because there is
    // no singularity in bessel_j for z>=0.
    return 0;
}

double minz_y(const int n)
{
    // return max(3., n)
    return .5;
}

double maxz_j(const int n)
{
    double z = (n * n + n + 1) / 1.221e-4;
    if (z >= 1000)
        z = std::max(1000, n * n);
    return z;
}

double maxz_y(const int n)
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

const table_type get_z_table(const double from, const double to, const double delta)
{
    table_type table;
    table.reserve(int((to - from) / delta));
    for (double z(from); z < to; z += delta)
        table.push_back(z);
    return table;
}

matrix zeros(const int length0, const int length1)
{
    matrix retval(length0, length1);
    for (int i(0); i < length0; i++)
        for (int j(0); j < length1; j++)
            retval(i, j) = 0;
    return retval;
}

void calculate_jns(const int n, const double z,
        double j[], double jdot[])
{
    double *jp = new double[n+2];
    gsl_sf_bessel_jl_array(n+1, z, jp);
    j[0] = jp[0];
    jdot[0] = -jp[1];
    for (int l(1); l < n+1; ++l)
    {
        j[l] = jp[l];
        jdot[l] = (l*jp[l-1] - (l+1)*jp[l+1])/(2*l+1);
    }
    delete[] jp;
}

void calculate_yns(const int n, const double z,
        double y[], double ydot[])
{
    double *yp = new double[n+2];
    gsl_sf_bessel_yl_array(n+1, z, yp);
    y[0] = yp[0];
    ydot[0] = -yp[1];
    for (int l(1); l < n+1; ++l)
    {
        y[l] = yp[l];
        ydot[l] = (n*yp[l-1] - (n+1)*yp[l])/(2*n+1);
    }
    delete[] yp;
}

void set_matrix(matrix& mat, const int i, const double array[], const int n)
{
    for (int index(0); index < n; ++index)
    {
        mat(i, index) = array[index];
    }
}

sjy_table jnyn(const int n, const int resolution)
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
        calculate_jns(n, z, j, jdot);
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
            calculate_yns(n, z, y, ydot);
        set_matrix(y_table, i, y, n+1);
        set_matrix(ydot_table, i, ydot, n+1);
        delete[] y;
        delete[] ydot;
    }

    j_table = trans(j_table);
    jdot_table = trans(jdot_table);
    y_table = trans(y_table);
    ydot_table = trans(ydot_table);

    return sjy_table(z_table, delta, j_table, jdot_table,
            y_table, ydot_table);
}

int searchsorted(const table_type& z_table, const double value)
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

valdot_type get_sub_sequence_from_matrix(const matrix mat0,
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

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Require output file name." << std::endl;
        return 1;
    }

    std::string filename(argv[1]);
    std::cout << "Generating as " << filename << std::endl;

    const int minn_j(4);
    const int maxn_j(51);
    const int minn_y(3);
    const int maxn_y(40);
    const int resolution(35);

    const std::string header_name("SPHERICAL_BESSEL_TABLE_HPP");
    const std::string ns_name("sb_table");
    std::ofstream ofs(filename.c_str());
    write_header(ofs, header_name, ns_name);

    const sjy_table table(jnyn(std::max(maxn_j, maxn_y), resolution));

    // j
    for (int n(minn_j); n<= maxn_j; ++n)
    {
        const int start(searchsorted(table.z, minz_j(n)));
        const double z_start(table.z.at(start));
        const int end(searchsorted(table.z, maxz_j(n)));
        const valdot_type js(get_sub_sequence_from_matrix(table.j, table.jdot, n, start, end));
        std::ostringstream oss0;
        oss0 << "sj_table" << n << "_f";
        write_arrays(ofs, oss0.str(), js);
        std::ostringstream oss1;
        oss1 << "sj_table" << n;
        write_table(ofs, oss1.str(), end-start, z_start, table.delta);
        if (n == maxn_j)
            std::cout << 'j' << js.size() << std::endl;
    }

    // y
    for (int n(minn_y); n<= maxn_y; ++n)
    {
        const int start(searchsorted(table.z, minz_y(n)));
        const double z_start(table.z.at(start));
        const int end(searchsorted(table.z, maxz_y(n)));
        const valdot_type ys(get_sub_sequence_from_matrix(table.y, table.ydot, n, start, end));
        std::ostringstream oss0;
        oss0 << "sy_table" << n << "_f";
        write_arrays(ofs, oss0.str(), ys);
        std::ostringstream oss1;
        oss1 << "sy_table" << n;
        write_table(ofs, oss1.str(), end-start, z_start, table.delta);
        if (n == maxn_y)
            std::cout << 'y' << ys.size() << std::endl;
    }

    write_table_array(ofs, "sj_table", minn_j, maxn_j);
    write_table_array(ofs, "sy_table", minn_y, maxn_y);

    write_footer(ofs, header_name, ns_name);

    ofs.close();

    return 0;
}
