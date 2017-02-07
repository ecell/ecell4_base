#ifndef ECELL_BARYCENTRIC
#define ECELL_BARYCENTRIC
#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Triangle.hpp>
#include <ostream>

namespace ecell4
{
namespace bd
{

template<typename realT>
struct Barycentric
{
    typedef realT       real_type;
    typedef real_type   value_type;
    typedef std::size_t size_type;
    typedef size_type   index_type;

    Barycentric(){}
    ~Barycentric(){}

    Barycentric(const real_type a, const real_type b, const real_type c);
    Barycentric(const Barycentric& b);
    Barycentric& operator=(const Barycentric& b);

    real_type  operator[](const index_type i) const {return val[i];}
    real_type& operator[](const index_type i)       {return val[i];}
    real_type  at(const index_type i) const {return val.at(i);}
    real_type& at(const index_type i)       {return val.at(i);}

  private:
    boost::array<realT, 3> val;
};

template<typename realT>
inline Barycentric<realT>::Barycentric(
        const real_type a, const real_type b, const real_type c)
{
    val[0] = a;
    val[1] = b;
    val[2] = c;
}

template<typename realT>
inline Barycentric<realT>::Barycentric(const Barycentric<realT>& b)
{
    val[0] = b.val[0];
    val[1] = b.val[1];
    val[2] = b.val[2];
}

template<typename realT>
inline Barycentric<realT>& Barycentric<realT>::operator=(
        const Barycentric<realT>& b)
{
    val[0] = b.val[0];
    val[1] = b.val[1];
    val[2] = b.val[2];
    return *this;
}

template<typename realT>
inline Barycentric<realT>
operator+(const Barycentric<realT>& lhs, const Barycentric<realT>& rhs)
{
    return Barycentric<realT>(lhs[0]+rhs[0], lhs[1]+rhs[1], lhs[2]+rhs[2]);
}


template<typename realT>
inline Barycentric<realT>
operator-(const Barycentric<realT>& lhs, const Barycentric<realT>& rhs)
{
    return Barycentric<realT>(lhs[0]-rhs[0], lhs[1]-rhs[1], lhs[2]-rhs[2]);
}

template<typename realT>
inline bool on_plane(
        const Barycentric<realT>& bary, const realT tolerance = 1e-10)
{
    return std::abs(bary[0] + bary[1] + bary[2] - 1.0) < tolerance;
}

template<typename realT>
inline bool is_inside(const Barycentric<realT>& bary)
{
    return on_plane(bary) &&
           (0. <= bary[0] && bary[0] <= 1.0) &&
           (0. <= bary[1] && bary[1] <= 1.0) &&
           (0. <= bary[2] && bary[2] <= 1.0);
}

template<typename realT>
inline Real3 to_absolute(const Barycentric<realT>& bary, const Triangle& tri)
{
    return tri.vertex_at(0) * bary[0] +
           tri.vertex_at(1) * bary[1] +
           tri.vertex_at(2) * bary[2];
}

namespace detail
{

template<typename realT>
inline realT triangle_area_2D(const realT x1, const realT y1,
        const realT x2, const realT y2, const realT x3, const realT y3)
{
    return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2);
}

} // detail

template<typename vectorT> //XXX normally, this is Real3.
Barycentric<typename vectorT::value_type>
to_barycentric(const vectorT& pos, const Triangle& face)
{
    typedef typename vectorT::value_type real_type;

    const vectorT&  a = face.vertex_at(0);
    const vectorT&  b = face.vertex_at(1);
    const vectorT&  c = face.vertex_at(2);
    const vectorT   m = cross_product(face.edge_at(0), face.edge_at(2)) * (-1.);
    const real_type x = std::abs(m[0]);
    const real_type y = std::abs(m[1]);
    const real_type z = std::abs(m[2]);

    real_type nu, nv, ood;
    if(x >= y && x >= z)
    {
        nu = detail::triangle_area_2D(pos[1], pos[2], b[1], b[2], c[1], c[2]);
        nv = detail::triangle_area_2D(pos[1], pos[2], c[1], c[2], a[1], a[2]);
        ood = 1.0 / m[0];
    }
    else if(y >= x && y >= z)
    {
        nu = detail::triangle_area_2D(pos[0], pos[2], b[0], b[2], c[0], c[2]);
        nv = detail::triangle_area_2D(pos[0], pos[2], c[0], c[2], a[0], a[2]);
        ood = 1.0 / -m[1];
    }
    else
    {
        nu = detail::triangle_area_2D(pos[0], pos[1], b[0], b[1], c[0], c[1]);
        nv = detail::triangle_area_2D(pos[0], pos[1], c[0], c[1], a[0], a[1]);
        ood = 1.0 / m[2];
    }
    Barycentric<real_type> bary;
    bary[0] = nu * ood;
    bary[1] = nv * ood;
    bary[2] = 1.0 - bary[0] - bary[1];
    return bary;
}

template<typename realT, typename charT, typename traitsT>
inline std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const Barycentric<realT>& b)
{
    os << b[0] << ", " << b[1] << ", " << b[2];
    return os;
}


} // bd
} // ecell4
#endif /* ECELL_BARYCENTRIC */
