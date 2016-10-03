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
#include "make_table_util.hpp"

#include "cjy_table.hpp"
#include "make_table_util.hpp"

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Require output file name." << std::endl;
        return 1;
    }

    std::string filename(argv[1]);
    std::cout << "Generating as " << filename << std::endl;

    const int minn_j(0);
    const int maxn_j(50);
    const int minn_y(0);
    const int maxn_y(50);
    const int resolution(35);

    const std::string header_name("CYLINDRICAL_BESSEL_TABLE_HPP");
    const std::string ns_name("cb_table");
    std::ofstream ofs(filename.c_str());
    write_header(ofs, header_name, ns_name);

    const cjy_table table(JnYn(std::max(maxn_j, maxn_y), resolution));

    // j
    for (int n(minn_j); n<= maxn_j; ++n)
    {
        const int start(searchsorted(table.z, minz_j(n)));
        const double z_start(table.z.at(start));
        const int end(searchsorted(table.z, maxz_j(n)));
        const valdot_type js(get_sub_sequence_from_matrix(table.j, table.jdot, n, start, end));
        std::ostringstream oss0;
        oss0 << "cj_table" << n << "_f";
        write_arrays(ofs, oss0.str(), js);
        std::ostringstream oss1;
        oss1 << "cj_table" << n;
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
        oss0 << "cy_table" << n << "_f";
        write_arrays(ofs, oss0.str(), ys);
        std::ostringstream oss1;
        oss1 << "cy_table" << n;
        write_table(ofs, oss1.str(), end-start, z_start, table.delta);
        if (n == maxn_y)
            std::cout << 'y' << ys.size() << std::endl;
    }

    write_table_array(ofs, "cj_table", minn_j, maxn_j);
    write_table_array(ofs, "cy_table", minn_y, maxn_y);

    write_footer(ofs, header_name, ns_name);

    ofs.close();

    return 0;
}
