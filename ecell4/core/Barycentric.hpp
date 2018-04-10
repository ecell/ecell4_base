#ifndef ECELL_BARYCENTRIC_HPP
#define ECELL_BARYCENTRIC_HPP
#include <ecell4/core/Triangle.hpp>

namespace ecell4
{

struct Barycentric
{
    typedef Real        value_type;
    typedef std::size_t size_type;
    typedef size_type   index_type;
    typedef boost::array<Real, 3> container_type;

    Barycentric(){}
    ~Barycentric(){}

    Barycentric(const Real a, const Real b, const Real c) throw();
    Barycentric(const Barycentric& b) throw();
    Barycentric& operator=(const Barycentric& b) throw();

    Real  operator[](const index_type i) const throw() {return val[i];}
    Real& operator[](const index_type i)       throw() {return val[i];}
    Real  at(const index_type i) const {return val.at(i);}
    Real& at(const index_type i)       {return val.at(i);}

  private:
    boost::array<Real, 3> val;
};

inline Barycentric::Barycentric(const Real a, const Real b, const Real c) throw()
{
    val[0] = a;
    val[1] = b;
    val[2] = c;
}

inline Barycentric::Barycentric(const Barycentric& b) throw()
{
    val[0] = b.val[0];
    val[1] = b.val[1];
    val[2] = b.val[2];
}

inline Barycentric& Barycentric::operator=(const Barycentric& b) throw()
{
    val[0] = b.val[0];
    val[1] = b.val[1];
    val[2] = b.val[2];
    return *this;
}

inline Barycentric
operator+(const Barycentric& lhs, const Barycentric& rhs) throw()
{
    return Barycentric(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
}

inline Barycentric
operator-(const Barycentric& lhs, const Barycentric& rhs) throw()
{
    return Barycentric(lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]);
}

inline bool on_plane(const Barycentric& bary, const Real tol = 1e-10) throw()
{
    return std::abs(bary[0] + bary[1] + bary[2] - 1.0) < tol;
}

inline bool is_inside(const Barycentric& bary) throw()
{
    return on_plane(bary) &&
           (0. <= bary[0] && bary[0] <= 1.0) &&
           (0. <= bary[1] && bary[1] <= 1.0) &&
           (0. <= bary[2] && bary[2] <= 1.0);
}

inline bool is_inside(const Barycentric& bary, const Real tolerance) throw()
{
    return on_plane(bary) &&
           (0.0 - tolerance <= bary[0] && bary[0] <= 1.0 + tolerance) &&
           (0.0 - tolerance <= bary[1] && bary[1] <= 1.0 + tolerance) &&
           (0.0 - tolerance <= bary[2] && bary[2] <= 1.0 + tolerance);
}

std::pair<std::size_t, Real>
first_cross_edge(const Barycentric& pos, const Barycentric& disp);

Barycentric force_put_inside(const Barycentric& bary);

Barycentric to_barycentric(const Real3& pos, const Triangle& face);

inline Real3 to_absolute(const Barycentric& bary, const Triangle& tri) throw()
{
    return tri.vertex_at(0) * bary[0] +
           tri.vertex_at(1) * bary[1] +
           tri.vertex_at(2) * bary[2];
}

template<typename charT, typename traitsT>
inline std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const Barycentric& b)
{
    os << b[0] << ", " << b[1] << ", " << b[2];
    return os;
}

}// ecell4
#endif // ECELL_BARYCENTRIC_HPP
