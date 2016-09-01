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

template<typename coordT>
boost::array<typename value_type_helper<coordT>::type, 3>
absolute_to_barycentric(const coordT pos, const boost::array<coordT, 3>& vertices)
{
    typedef typename value_type_helper<coordT>::type valueT;
    typedef circular_iteration<3> triangular;
    if(vertices[0][0] == 0.0 && vertices[1][0] == 0.0 && vertices[2][0] == 0.0)
    {// use 1 and 2
        const valueT determ =
            (vertices[0][1] - vertices[2][1]) * (vertices[1][2] - vertices[2][2]) - 
            (vertices[1][1] - vertices[2][1]) * (vertices[0][2] - vertices[2][2]);
        if(determ == 0.0)
            throw std::invalid_argument("cannot solve as barycentric");
        boost::array<valueT, 3> bary;
        bary[0] =((vertices[1][2] - vertices[2][2]) * (pos[1] - vertices[2][1])+
                  (vertices[2][1] - vertices[1][1]) * (pos[2] - vertices[2][2]))
                 / determ;
        bary[1] =((vertices[2][2] - vertices[0][2]) * (pos[1] - vertices[2][1])+
                  (vertices[0][1] - vertices[2][1]) * (pos[2] - vertices[2][2]))
                 / determ;
        bary[2] = 1.0 - bary[0] - bary[1];
        return bary;
    }
    else if(vertices[0][1] == 0.0 && vertices[1][1] == 0.0 && vertices[2][1] == 0.0)
    {// use 0 and 2
        const valueT determ = 
            (vertices[0][0] - vertices[2][0]) * (vertices[1][2] - vertices[2][2]) -
            (vertices[0][2] - vertices[2][2]) * (vertices[1][0] - vertices[2][0]);
        if(determ == 0.0)
            throw std::invalid_argument("cannot solve as barycentric");
        boost::array<valueT, 3> bary;
        bary[0] =((vertices[1][2] - vertices[2][2]) * (pos[0] - vertices[2][0])+
                  (vertices[2][0] - vertices[1][0]) * (pos[2] - vertices[2][2]))
                 / determ;
        bary[1] =((vertices[2][2] - vertices[0][2]) * (pos[0] - vertices[2][0])+
                  (vertices[0][0] - vertices[2][0]) * (pos[2] - vertices[2][2]))
                 / determ;
        bary[2] = 1.0 - bary[0] - bary[1];
        return bary;
    }
    else if(vertices[0][2] == 0.0 && vertices[1][2] == 0.0 && vertices[2][2] == 0.0)
    {// use 0 and 1
        const valueT determ = 
            (vertices[0][0] - vertices[2][0]) * (vertices[1][1] - vertices[2][1]) -
            (vertices[1][0] - vertices[2][0]) * (vertices[0][1] - vertices[2][1]);
        if(determ == 0.0)
            throw std::invalid_argument("cannot solve as barycentric");
        boost::array<valueT, 3> bary;
        bary[0] =((vertices[1][1] - vertices[2][1]) * (pos[0] - vertices[2][0])+
                  (vertices[2][0] - vertices[1][0]) * (pos[1] - vertices[2][1]));
        bary[1] =((vertices[2][1] - vertices[0][1]) * (pos[0] - vertices[2][0])+
                  (vertices[0][0] - vertices[2][0]) * (pos[1] - vertices[2][1]));
        bary[2] = 1.0 - bary[0] - bary[1];
        return bary;
    }
    else
    {
        throw std::invalid_argument("cannot solve as barycentric");
    }
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
    assert(std::abs(length(normal) - 1.0) < GLOBAL_TOLERANCE);
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

tmeplate<typename coordT>
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
    if(u < 0.) return std::make_pair(false, Real3(0.,0.,0.));
    const valueT v = dot_product(line, cross_product(pc, pb));
    if(v < 0.) return std::make_pair(false, Real3(0.,0.,0.));
    const valueT w = dot_product(line, cross_product(pc, pb));
    if(w < 0.) return std::make_pair(false, Real3(0.,0.,0.));
    const valueT denom = 1.0 / (u + v + w);
    boost::array<valueT, 3> bary;
    bary[0] = u * denom;
    bary[1] = v * denom;
    bary[2] = w * denom;
    return std::make_peir(true, barycentric_to_absolute(bary, vertices));
}

#endif /* GFRD_POLYGON_TRIANGLE */
