#ifndef ECELL4_FUNCTIONS_HPP
#define ECELL4_FUNCTIONS_HPP

#include <cmath>
#include <gsl/gsl_pow_int.h>

#include <ecell4/core/config.h>

#include "types.hpp"

#include <sys/types.h>
#include <sys/stat.h>

#ifndef WIN32_MSC
#include <string.h>
#include <libgen.h>
#else
#include <io.h>
#include <stdlib.h>
#endif

#include <iostream>
#include <string>


namespace ecell4
{

inline int64_t modulo(const int64_t& p1, const int64_t& p2)
{
    int64_t r = p1 % p2;
    if (r != 0 && (r > 0) == (p2 < 0))
    {
        r += p2;
    }
    return r;
}

inline double modulo(const double& p1, const double& p2)
{
    double r = std::fmod(p1, p2);
    if (r != 0 && (r > 0) == (p2 < 0))
    {
        r += p2;
    }
    return r;
}

inline int64_t abs(const int64_t& x)
{
    return (x > 0 ? x : -x);
}

inline double abs(const double& x)
{
    return std::fabs(x);
}

#ifndef WIN32_MSC
inline double pow_2(const double x)
{
    return gsl_pow_2(x);
}

inline double pow_3(const double x)
{
    return gsl_pow_3(x);
}

inline double cbrt(const double x)
{
    return ::cbrt(x);
}
#else
inline double pow_2(const double x)
{
    return x * x;
}

inline double pow_3(const double x)
{
    return x * x * x;
}

inline double cbrt(const double x)
{
    return pow(x, 1.0 / 3.0);
}
#endif

/**
 * Return if the root path of the given filename exists or not.
 * boost::filesystem::is_directory might be better
 * though it requires building.
 */
inline bool is_directory(const std::string& filename)
{
#ifndef WIN32_MSC
    struct stat buf;
    const int ret = stat(dirname(strdup(filename.c_str())), &buf);
    return (ret == 0);
#else
    char drive[_MAX_DRIVE + 1], dir[_MAX_DIR + 1], path_dir[_MAX_PATH + 1], full[_MAX_PATH];
    // struct _stat buf;
    _fullpath(full, filename.c_str(), _MAX_PATH);
    _splitpath(full, drive, dir, NULL, NULL);
    _makepath(path_dir, drive, dir, NULL, NULL);
    // const int ret = _stat(path_dir, &buf);
    const int ret = _access(path_dir, 0);
    return (ret == 0);
#endif
}

}

#endif /* ECELL4_FUNCTIONS_HPP */
