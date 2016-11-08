#ifndef GFRD_POLYGON_BARYCENTRIC
#define GFRD_POLYGON_BARYCENTRIC
#include <boost/array.hpp>

template<typename realT>
struct Barycentric
{
    typedef realT real_type;
    typedef real_type value_type;

    boost::array<value_type, 3> coordinate;

    Barycentric(){}
    ~Barycentric(){}

    explicit Barycentric(const boost::array<value_type, 3>& coord)
        : coordinate(coord){}

    Barycentric(const value_type a, const value_type b, const value_type c)
    {
        coordinate[0] = a;
        coordinate[1] = b;
        coordinate[2] = c;
    }
};

template<typename realT>
inline realT
triangle_area_2D(const realT x1, const realT y1, const realT x2, const realT y2,
                 const realT x3, const realT y3)
{
    return (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2);
}

template<typename coordT>
Barycentric<typename element_type_of<coordT>::type>
make_barycentric(const coordT& pos, const boost::array<coordT, 3>& tri)
{
    // the implementation of this function is based on
    // Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // p. 51

    typedef typename element_type_of<coordT>::type scalarT;
    const coordT& a = tri.at(0);
    const coordT& b = tri.at(1);
    const coordT& c = tri.at(2);
    const coordT m = cross_product(b - a, c - a);
    const scalarT x = std::abs(m[0]);
    const scalarT y = std::abs(m[1]);
    const scalarT z = std::abs(m[2]);

    scalarT nu, nv, ood;
    if (x >= y && x >= z)
    {
        nu = triangle_area_2D(pos[1], pos[2], b[1], b[2], c[1], c[2]);
        nv = triangle_area_2D(pos[1], pos[2], c[1], c[2], a[1], a[2]);
        ood = 1.0 / m[0];
    }
    else if (y >= x && y >= z)
    {
        nu = triangle_area_2D(pos[0], pos[2], b[0], b[2], c[0], c[2]);
        nv = triangle_area_2D(pos[0], pos[2], c[0], c[2], a[0], a[2]);
        ood = -1.0 / m[1];
    }
    else
    {
        nu = triangle_area_2D(pos[0], pos[1], b[0], b[1], c[0], c[1]);
        nv = triangle_area_2D(pos[0], pos[1], c[0], c[1], a[0], a[1]);
        ood = 1.0 / m[2];
    }
    Barycentric<scalarT> bary;
    bary.coordinate[0] = nu * ood;
    bary.coordinate[1] = nv * ood;
    bary.coordinate[2] = 1.0 - bary.coordinate[0] - bary.coordinate[1];
    return bary;
}

template<typename coordT>
inline coordT make_absolute(
    const Barycentric<typename element_type_of<coordT>::type>& barycentric,
    const boost::array<coordT, 3>& triangle)
{
    return triangle.at(0) * barycentric.coordinate[0] +
           triangle.at(1) * barycentric.coordinate[1] +
           triangle.at(2) * barycentric.coordinate[2];
}

template<typename realT>
inline bool is_inside(const Barycentric<realT>& barycentric)
{
    return
        (0. <= barycentric.coordinate[0] && barycentric.coordinate[0] <= 1.0) &&
        (0. <= barycentric.coordinate[1] && barycentric.coordinate[1] <= 1.0) &&
        (0. <= barycentric.coordinate[2] && barycentric.coordinate[2] <= 1.0);
}



#endif /* SURFER_BARYCENTRIC */
