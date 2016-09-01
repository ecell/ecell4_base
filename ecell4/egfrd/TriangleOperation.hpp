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
distance(const coordT& pos,
    const boost::array<coordT, 3>& vertices,
    const boost::array<coordT, 3>& edges,
    const boost::array<typename value_type_helper<coordT>::type, 3>& side_length,
    const coordT& normal)
{
    typedef typename value_type_helper<coordT>::type valueT;
    const coordT projected = project_to_plane(pos, vertices, normal);
    const boost::array<valueT, 3> barycentric = 
        absolute_to_barycentric(pos, vertices);

    // never be (-, -, -) because projected pos and triangle are on same plane
    if(barycentric[0] < 0 && barycentric[1] < 0 && barycentric[2] < 0)
        throw std::logic_error("projection error");

    if(barycentric[0] >= 0 && barycentric[1] >= 0 && barycentric[2] >= 0)
    {// nearest position is inside of triangle. r != 0
        const valueT s =
            (side_length[0] + side_length[1] + side_length[2]) * 0.5;
        const valueT area = std::sqrt(s *
                (s - side_length[0]) * (s - side_length[1]) * (s - side_length[2]));
        boost::array<valueT, 3> perpendicular;
        perpendicular[0] = 2.0 * area * barycentric[0] / side_length[1];
        perpendicular[1] = 2.0 * area * barycentric[1] / side_length[2];
        perpendicular[2] = 2.0 * area * barycentric[2] / side_length[0];

        return std::make_pair(
                length(pos - barycentric_to_absolute(barycentric, vertices)),
                (*std::min_element(barycentric.begin(), barycentric.end())));
    }
    else if(barycentric[0] * barycentric[1] * barycentric[2] > 0) //(+, -, -)
    {// nearest position is vertex
        if(barycentric[0] > 0)
            return std::make_pair(length(pos - vertices[0]), 0);
        else if(barycentric[1] > 0)
            return std::make_pair(length(pos - vertices[1]), 0);
        else if(barycentric[2] > 0)
            return std::make_pair(length(pos - vertices[2]), 0);
        else
            throw std::logic_error(
                    "distance between point to triangle: never reach here");
    }
    else // (+, -, -)
    {//nearest position is on edge, or vertex
        if(barycentric[0] < 0.0)
        {
            const valueT tmp =
                dot_product(edges[1], (pos - vertices[1])) / length_sq(edges[1]);
            if(0.0 < tmp && tmp < 1.0)
            {
                return std::make_pair(
                        length(pos - (vertices[1] + edges[1] * tmp)), 0.0);
            }
            else if(tmp <= 0.0)
            {
                return std::make_pair(length(pos - vertices[1]), 0);
            }
            else if(1.0 <= tmp)
            {
                return std::make_pair(length(pos - vertices[2]), 0);
            }
        }
        else if(barycentric[1] < 0)
        {
            const valueT tmp =
                dot_product(edges[2], (pos - vertices[2])) / length_sq(edges[2]);
            if(0.0 < tmp && tmp < 1.0)
            {
                return std::make_pair(
                        length(pos - (vertices[2] + edges[2] * tmp)), 0.0);
            }
            else if(tmp <= 0.0)
            {
                return std::make_pair(length(pos - vertices[2]), 0);
            }
            else if(1.0 <= tmp)
            {
                return std::make_pair(length(pos - vertices[0]), 0);
            }
        }
        else if(barycentric[2] < 0)
        {
            const valueT tmp =
                dot_product(edges[0], (pos - vertices[0])) / length_sq(edges[0]);
            if(0.0 < tmp && tmp < 1.0)
            {
                return std::make_pair(
                        length(pos - (vertices[0] + edges[0] * tmp)), 0.0);
            }
            else if(tmp <= 0.0)
            {
                return std::make_pair(length(pos - vertices[0]), 0);
            }
            else if(1.0 <= tmp)
            {
                return std::make_pair(length(pos - vertices[1]), 0);
            }
        }
        else
        {
            throw std::logic_error(
                    "distance between point to triangle: never reach here");
        }
    }
    throw std::logic_error(
            "distance between point to triangle: never reach here");
}

template<typename coordT>
std::pair<bool, coordT> // pair of (whether pierce), pierce point
is_pierce(const coordT& begin, const coordT& end,
    const boost::array<coordT, 3>& vertices)
{
    typedef typename value_type_helper<coordT>::type valueT;
    const coordT line = end - begin;
    const coordT pa   = vertices[0] - begin;
    const coordT pb   = vertices[1] - begin;
    const coordT pc   = vertices[2] - begin;
    const valueT u = dot_product(line, cross_product(pc, pb));
    if(u < 0.) return std::make_pair(false, coordT(0.,0.,0.));
    const valueT v = dot_product(line, cross_product(pa, pc));
    if(v < 0.) return std::make_pair(false, coordT(0.,0.,0.));
    const valueT w = dot_product(line, cross_product(pb, pa));
    if(w < 0.) return std::make_pair(false, coordT(0.,0.,0.));
    const valueT denom = 1.0 / (u + v + w);
    boost::array<valueT, 3> bary;
    bary[0] = u * denom;
    bary[1] = v * denom;
    bary[2] = w * denom;
    return std::make_pair(true, barycentric_to_absolute(bary, vertices));
}

#endif /* GFRD_POLYGON_TRIANGLE */
