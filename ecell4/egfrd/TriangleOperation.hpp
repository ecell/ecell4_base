#ifndef GFRD_POLYGON_TRIANGLE_OPERATION
#define GFRD_POLYGON_TRIANGLE_OPERATION
#include "Real3Type.hpp"
#include <boost/array.hpp>
#include <algorithm>
#include <cassert>

namespace ecell4{

template<typename coordT>
inline coordT centroid(const boost::array<coordT, 3>& vertices)
{
    return (vertices[0] + vertices[1] + vertices[2]) * (1e0 / 3e0);
}

template<typename coordT>
inline coordT incenter(const boost::array<coordT, 3>& vertices)
{
    typedef typename element_type_of<coordT>::type valueT;
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
    typedef typename element_type_of<coordT>::type valueT;
    const valueT a = length(edges[1]);
    const valueT b = length(edges[2]);
    const valueT c = length(edges[0]);
    const valueT abc = a + b + c;
    return (vertices[0] * a + vertices[1] * b + vertices[2] * c) * (1e0 / abc);
}

template<typename coordT>
inline coordT incenter(const boost::array<coordT, 3>& vertices,
    const boost::array<typename element_type_of<coordT>::type, 3>& length_of_edge)
{
    typedef typename element_type_of<coordT>::type valueT;
    const valueT a = length_of_edge[1];
    const valueT b = length_of_edge[2];
    const valueT c = length_of_edge[0];
    const valueT abc = a + b + c;
    return (vertices[0] * a + vertices[1] * b + vertices[2] * c) * (1e0 / abc);
}

template<typename coordT>
inline std::size_t
match_edge(const coordT& vec, const boost::array<coordT, 3>& edges,
           const typename element_type_of<coordT>::type tol = 1e-10)
{
    for(std::size_t i=0; i<3; ++i)
    {
        if((std::abs(vec[0] - edges[i][0]) < tol) &&
           (std::abs(vec[1] - edges[i][1]) < tol) &&
           (std::abs(vec[2] - edges[i][2]) < tol)) return i;
    }
    throw std::invalid_argument("not match any edge");
}

template<typename coordT>
coordT
project_to_plane(const coordT& pos, const boost::array<coordT, 3>& vertices,
                 const coordT& normal)
{
    typedef typename element_type_of<coordT>::type valueT;
    assert(std::abs(length(normal) - 1.0) < 1e-10);
    const valueT distance = dot_product(normal, pos - vertices.front());
    return pos - (normal * distance);
}

template<typename coordT>
coordT
closest_point(const coordT& pos, const boost::array<coordT, 3>& vertices)
{
    typedef typename element_type_of<coordT>::type valueT;
    // this implementation is from Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // pp.141-142

    const coordT a = vertices[0];
    const coordT b = vertices[1];
    const coordT c = vertices[2];

    const coordT ab = b - a;
    const coordT ac = c - a;
    const coordT ap = pos - a;
    const valueT d1 = dot_product(ab, ap);
    const valueT d2 = dot_product(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0)
        return a;

    const coordT bp = pos - b;
    const valueT d3 = dot_product(ab, bp);
    const valueT d4 = dot_product(ac, bp);
    if (d3 >= 0.0 && d4 <= d3)
        return b;

    const valueT vc = d1*d4 - d3*d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        valueT v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const coordT cp = pos - c;
    const valueT d5 = dot_product(ab, cp);
    const valueT d6 = dot_product(ac, cp);
    if (d6 >= 0.0 && d5 <= d6)
        return c;

    const valueT vb = d5*d2 - d1*d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const valueT w = d2 / (d2 - d6);
        return a + ac * w;
    }

    const valueT va = d3*d6 - d5*d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const valueT w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }

    const valueT denom = 1.0 / (va + vb + vc);
    const valueT v = vb * denom;
    const valueT w = vc * denom;
    return a + ab * v + ac * w;
}

template<typename coordT>
std::pair<typename element_type_of<coordT>::type, // distance
          typename element_type_of<coordT>::type> // r of circle in triangle
distance(const coordT& pos, const boost::array<coordT, 3>& vertices)
{
    return std::make_pair(length(closest_point(pos, vertices) - pos), 0.);
}

template<typename coordT>
std::pair<bool, coordT>
test_intersect_segment_triangle(const coordT& begin, const coordT& end,
                                const boost::array<coordT, 3>& vertices)
{
    typedef typename element_type_of<coordT>::type valueT;
    // this implementation is from Real-Time Collision Detection by Christer Ericson,
    // published by Morgan Kaufmann Publishers, (c) 2005 Elsevier Inc.
    // pp.190-194

    const coordT line = begin - end;
    const coordT ab = vertices[1] - vertices[0];
    const coordT ac = vertices[2] - vertices[0];
    const coordT normal = cross_product(ab, ac);

    const valueT d = dot_product(line, normal);
    if(d < 0.0)
        return std::make_pair(false, coordT(0.,0.,0.));

    const coordT ap = begin - vertices[0];
    const valueT t = dot_product(ap, normal);
    if(t < 0.0 || d < t)
        return std::make_pair(false, coordT(0.,0.,0.));

    const coordT e = cross_product(line, ap);
    valueT v = dot_product(ac, e);
    if(v < 0. || d < v)
        return std::make_pair(false, coordT(0.,0.,0.));
    valueT w = -1.0 * dot_product(ab, e);
    if(w < 0. || d < v + w)
        return std::make_pair(false, coordT(0.,0.,0.));

    const valueT ood = 1. / d;
    v *= ood;
    w *= ood;
    const valueT u = 1. - v - w;
    const coordT intersect = vertices[0] * u + vertices[1] * v + vertices[2] * w;

    return std::make_pair(true, intersect);
}


}
#endif /* GFRD_POLYGON_TRIANGLE */
