#ifndef GFRD_POLYGON_TRIANGLE_OPERATION
#define GFRD_POLYGON_TRIANGLE_OPERATION
#include "Vector3Operation.hpp"
#include "circular_iteration.hpp"
#include <boost/array.hpp>
#include <algorithm>
#include <cassert>

template<typename coordT>
inline coordT centroid(const boost::array<coordT, 3>& vertices)
{
    return (vertices[0] + vertices[1] + vertices[2]) * (1e0 / 3e0);
}

template<typename coordT>
inline coordT incenter(const boost::array<coordT, 3>& vertices)
{
    typedef typename value_type_helper<coordT>::type valueT;
    const valueT a = length(vertices[2] - vertices[1]);
    const valueT b = length(vertices[0] - vertices[2]);
    const valueT c = length(vertices[1] - vertices[0]);
    const valueT abc = a + b + c;
    return (vertices[0] * a + vertices[1] * b + vertices[2] * c) * (1e0 / abc);
}

template<typename coordT>
inline coordT incenter(const boost::array<coordT, 3>& vertices,
                       const boost::array<coordT, 3>& edges)
{
    typedef typename value_type_helper<coordT>::type valueT;
    const valueT a = length(edges[1]);
    const valueT b = length(edges[2]);
    const valueT c = length(edges[0]);
    const valueT abc = a + b + c;
    return (vertices[0] * a + vertices[1] * b + vertices[2] * c) * (1e0 / abc);
}

template<typename coordT>
inline coordT incenter(const boost::array<coordT, 3>& vertices,
    const boost::array<typename value_type_helper<coordT>::type, 3>& length_of_edge)
{
    typedef typename value_type_helper<coordT>::type valueT;
    const valueT a = length_of_edge[1];
    const valueT b = length_of_edge[2];
    const valueT c = length_of_edge[0];
    const valueT abc = a + b + c;
    return (vertices[0] * a + vertices[1] * b + vertices[2] * c) * (1e0 / abc);
}

template<typename coordT>
inline std::size_t match_edge(const coordT& vec,
                              const boost::array<coordT, 3>& edges)
{
    for(std::size_t i=0; i<3; ++i)
        if(is_same_vec(vec, edges[i])) return i;
    throw std::invalid_argument("not match any edge");
}

template<typename realT>
inline realT
triangle_area_2D(const realT x1, const realT y1, const realT x2, const realT y2,
                 const realT x3, const realT y3)
{
    return (x1-x2)*(y2-y3) - (x2-x3)*(y1-y2);
}

template<typename coordT>
boost::array<typename value_type_helper<coordT>::type, 3>
absolute_to_barycentric(const coordT pos, const boost::array<coordT, 3>& vertices)
{
    // this implementation is from Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.

    typedef typename value_type_helper<coordT>::type valueT;
    const coordT& a = vertices[0];
    const coordT& b = vertices[1];
    const coordT& c = vertices[2];
    const coordT m = cross_product(b - a, c - a);
    // Nominators and one-over-denominator for u and v ratios
    valueT nu, nv, ood;
    // Absolute components for determining projection plane
    const valueT x = std::abs(m[0]);
    const valueT y = std::abs(m[1]);
    const valueT z = std::abs(m[2]);
    // Compute areas in plane of largest projection
    if (x >= y && x >= z) // x is largest, project to the yz plane
    {
        nu = triangle_area_2D(pos[1], pos[2], b[1], b[2], c[1], c[2]); // Area of PBC in yz plane
        nv = triangle_area_2D(pos[1], pos[2], c[1], c[2], a[1], a[2]); // Area of PCA in yz plane
        ood = 1.0 / m[0]; // 1/(2*area of ABC in yz plane)
    }
    else if (y >= x && y >= z)// y is largest, project to the xz plane
    {
        nu = triangle_area_2D(pos[0], pos[2], b[0], b[2], c[0], c[2]);
        nv = triangle_area_2D(pos[0], pos[2], c[0], c[2], a[0], a[2]);
        ood = 1.0 / -m[1];
    }
    else // z is largest, project to the xy plane
    {
        nu = triangle_area_2D(pos[0], pos[1], b[0], b[1], c[0], c[1]);
        nv = triangle_area_2D(pos[0], pos[1], c[0], c[1], a[0], a[1]);
        ood = 1.0 / m[2];
    }
    boost::array<valueT, 3> bary;
    bary[0] = nu * ood;
    bary[1] = nv * ood;
    bary[2] = 1.0 - bary[0] - bary[1];
    return bary;
}

template<typename coordT>
coordT barycentric_to_absolute(
        const boost::array<typename value_type_helper<coordT>::type,3>& bary,
        const boost::array<coordT, 3>& vertices)
{
    return vertices[0] * bary[0] + vertices[1] * bary[1] + vertices[2] * bary[2];
}

template<typename coordT>
coordT
project_to_plane(const coordT& pos, const boost::array<coordT, 3>& vertices,
                 const coordT& normal)
{
    typedef typename value_type_helper<coordT>::type valueT;
    assert(std::abs(length(normal) - 1.0) < 1e-10);
    const valueT distance = dot_product(normal, pos - vertices.front());
    return pos - (normal * distance);
}

template<typename coordT>
std::pair<typename value_type_helper<coordT>::type, // distance
          typename value_type_helper<coordT>::type> // r of circle in triangle
distance(const coordT& pos, const boost::array<coordT, 3>& vertices)
{
    typedef typename value_type_helper<coordT>::type valueT;
    // this implementation is from Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.

    coordT const& a = vertices[0];
    coordT const& b = vertices[1];
    coordT const& c = vertices[2];

    // Check if P in vertex region outside A
    coordT ab = b - a;
    coordT ac = c - a;
    coordT ap = pos - a;
    valueT d1 = dot_product(ab, ap);
    valueT d2 = dot_product(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0)
        return std::make_pair(length(pos - a), 0.0); // barycentric coordinates (1,0,0)
    // Check if P in vertex region outside B
    coordT bp = pos - b;
    valueT d3 = dot_product(ab, bp);
    valueT d4 = dot_product(ac, bp);
    if (d3 >= 0.0 && d4 <= d3)
        return std::make_pair(length(pos - b), 0.0); // barycentric coordinates (0,1,0)
    // Check if P in edge region of AB, if so return projection of P onto AB
    valueT vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        valueT v = d1 / (d1 - d3);
        return std::make_pair(length(a + ab * v- pos), 0.0); // barycentric coordinates (1-v,v,0)
    }
    // Check if P in vertex region outside C
    coordT cp = pos - c;
    valueT d5 = dot_product(ab, cp);
    valueT d6 = dot_product(ac, cp);
    if (d6 >= 0.0 && d5 <= d6)
        return std::make_pair(length(c - pos), 0.0); // barycentric coordinates (0,0,1)

    //Check if P in edge region of AC, if so return projection of P onto AC
    valueT vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        valueT w = d2 / (d2 - d6);
        return std::make_pair(length(a + ac * w - pos), 0.0); // barycentric coordinates (1-w,0,w)
    }
    // Check if P in edge region of BC, if so return projection of P onto BC
    valueT va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        valueT w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return std::make_pair(length(b + (c - b) * w - pos), 0.0); // barycentric coordinates (0,1-w,w)
    }
    // P inside face region. Compute Q through its barycentric coordinates (u,v,w)
    valueT denom = 1.0f / (va + vb + vc);
    valueT v = vb * denom;
    valueT w = vc * denom;
    return std::make_pair(length(a + ab * v + ac * w - pos), 0.0); // = u*a + v*b + w*c, u = va * denom = 1.0f-v-w
}

template<typename coordT>
std::pair<bool, coordT> // pair of (whether pierce), pierce point
is_pierce(const coordT& begin, const coordT& end,
          const boost::array<coordT, 3>& vertices)
{
    typedef typename value_type_helper<coordT>::type valueT;
    // this implementation is from Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    const coordT line = end - begin;
//     std::cerr << "line = " << line << std::endl;
    const coordT pa   = vertices[0] - begin;
//     std::cerr << "pa = " << pa << std::endl;
    const coordT pb   = vertices[1] - begin;
//     std::cerr << "pb = " << pb << std::endl;
    const coordT pc   = vertices[2] - begin;
//     std::cerr << "pc = " << pc << std::endl;
    const valueT u = dot_product(line, cross_product(pc, pb));
//     std::cerr << u << std::endl;
    if(u < 0.) return std::make_pair(false, coordT(0.,0.,0.));
    const valueT v = dot_product(line, cross_product(pa, pc));
//     std::cerr << v << std::endl;
    if(v < 0.) return std::make_pair(false, coordT(0.,0.,0.));
    const valueT w = dot_product(line, cross_product(pb, pa));
//     std::cerr << w << std::endl;
    if(w < 0.) return std::make_pair(false, coordT(0.,0.,0.));
    const valueT denom = 1.0 / (u + v + w);
    boost::array<valueT, 3> bary;
    bary[0] = u * denom;
    bary[1] = v * denom;
    bary[2] = w * denom;
    const coordT intersect = barycentric_to_absolute(bary, vertices);
//     std::cerr << "int = " << intersect << std::endl;
//     std::cerr << "intersect - begin = " << intersect - begin << std::endl;
    const valueT len_l = length(line);
    const valueT len_p = length(intersect - begin);
//     std::cerr << "len_l " << len_l << std::endl;
//     std::cerr << "len_p " << len_p << std::endl;
    const bool is_intersect = (len_l > len_p) && (dot_product(line, intersect - begin) > 0);

    return std::make_pair(is_intersect, intersect);
}

#endif /* GFRD_POLYGON_TRIANGLE */
