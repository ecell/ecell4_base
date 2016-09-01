#ifndef GFRD_POLYGON_PARAMETRIC_POSITION
#define GFRD_POLYGON_PARAMETRIC_POSITION
#include "Vector3Operation.hpp"
#include "circular_iteration.hpp"
#include <boost/array.hpp>
#include <algorithm>
#include <iostream>

// parametric representation of position in triangle
// this position is always relative.
/*      /\
 *  ^  /  \
 * b| /_.  \
 *   /_/____\
 *  o -->  
 *     a    */
template<typename T_value>
struct ParametricPosition
{
    typedef T_value value_type;
    ParametricPosition(){}
    ParametricPosition(const value_type a_, const value_type b_): a(a_), b(b_){}
    value_type a;
    value_type b;
};

template<typename T>
inline bool is_in_face(const ParametricPosition<T>& pos,
                       const T t = 1e-12)
{
    return (0e0 - t <= pos.a         && pos.a         <= 1e0 + t) &&
           (0e0 - t <=         pos.b &&         pos.b <= 1e0 + t) &&
           (0e0 - t <= pos.a + pos.b && pos.a + pos.b <= 1e0 + t);
}

template<typename T>
inline bool is_on_edge(const ParametricPosition<T>& pos,
                       const T t = 1e-12)
{
    return (0e0 - t <= pos.a         && pos.a         <= 0e0 + t) ||
           (0e0 - t <=         pos.b &&         pos.b <= 0e0 + t) ||
           (1e0 - t <= pos.a + pos.b && pos.a + pos.b <= 1e0 + t);
}

// template<typename T, typename T_edgeid = std::size_t>
template<typename T>
inline std::size_t which_edge(const ParametricPosition<T>& pos,
                              const T t = 1e-12)
{
         if(0e0 - t <= pos.a         && pos.a         <= 0e0 + t) return 2;
    else if(0e0 - t <=         pos.b &&         pos.b <= 0e0 + t) return 0;
    else if(1e0 - t <= pos.a + pos.b && pos.a + pos.b <= 1e0 + t) return 1;
    else throw std::invalid_argument("not on edge");
}

template<typename T>
inline bool is_on_vertex(const ParametricPosition<T>& pos, const T t = 0e0)
{
    return ((0e0-t <= pos.a && pos.a <= 0e0+t) && (0e0-t <= pos.b && pos.b <= 0e0+t)) ||
           ((0e0-t <= pos.a && pos.a <= 0e0+t) && (1e0-t <= pos.b && pos.b <= 1e0+t)) ||
           ((0e0-t <= pos.b && pos.b <= 0e0+t) && (1e0-t <= pos.a && pos.a <= 1e0+t));
}

template<typename T>
inline std::size_t which_vertex(const ParametricPosition<T>& pos, const T t = 0e0)
{
         if((0e0 - t <= pos.a && pos.a <= 0e0 + t) &&
            (0e0 - t <= pos.b && pos.b <= 0e0 + t)) return 0;
    else if((1e0 - t <= pos.a && pos.a <= 1e0 + t) &&
            (0e0 - t <= pos.b && pos.b <= 0e0 + t)) return 1;
    else if((0e0 - t <= pos.a && pos.a <= 0e0 + t) &&
            (1e0 - t <= pos.b && pos.b <= 1e0 + t)) return 2;
    else throw std::invalid_argument("not on vertex");
}

template<typename T>
T cross_ratio(const ParametricPosition<T>& pos,
              const ParametricPosition<T>& dis, const std::size_t edge_id)
{
    switch(edge_id)
    {
        case 0:
            if(dis.b == 0e0)         throw std::logic_error("divide by zero");
            return (-pos.b / dis.b);
        case 1:
            if(dis.a + dis.b == 0e0) throw std::logic_error("divide by zero");
            return (1e0 - pos.a - pos.b) / (dis.a + dis.b);
        case 2:
            if(dis.a == 0e0)         throw std::logic_error("divide by zero");
            return -pos.a / dis.a;
        default:
            throw std::invalid_argument("cross_ratio: invalid edge_id");
    }
}

template<typename T>
std::size_t cross_edge_lhs_is_in_face_rhs_is_out_of_face(
        const ParametricPosition<T>& lhs, const ParametricPosition<T>& rhs)
{
    if(is_on_edge(lhs, 0e0))
    {
        switch(which_edge(lhs, 0e0))
        {
            case 0:
                if(rhs.b < 0e0) return 0;
            case 1:
                if(rhs.a + rhs.b > 1e0) return 1;
            case 2:
                if(rhs.b < 0e0) return 2;
            default:
                throw std::logic_error("invalid local edge id");
        }
    }
/* case    */
/*   \5/   */
/*    v    */
/*  2 ^ 1  */
/* __/_\__ */
/* 3/ 0 \4 */
/* /     \ */
    // case 0, 1, or 2
    if((0e0 <= rhs.a && rhs.b < 0e0) && (rhs.a + rhs.b <= 1e0))
        return 0;
    else if((0e0 <= rhs.a && 0e0 <= rhs.b) && (1e0 < rhs.a + rhs.b))
        return 1;
    else if((rhs.a < 0e0 && 0e0 <= rhs.b) && (rhs.a + rhs.b <= 1e0))
        return 2;

    const T da = rhs.a - lhs.a;
    const T db = rhs.b - lhs.b;
    const ParametricPosition<T> disp(da, db);
    if(rhs.a < 0e0 && rhs.b < 0e0) // case 3
    {
        const T ratio0 = cross_ratio(lhs, disp, 0);
        const T ratio2 = cross_ratio(lhs, disp, 2);
             if(ratio0 < ratio2) return 0;
        else if(ratio0 > ratio2) return 2;
        else throw std::logic_error("displacement through just on vertex");
    }
    else if((rhs.a > 0e0 && rhs.b < 0e0) && (rhs.a + rhs.b > 1e0)) // case 4
    {
        const T ratio0 = cross_ratio(lhs, disp, 0);
        const T ratio1 = cross_ratio(lhs, disp, 1);
             if(ratio0 < ratio1) return 0;
        else if(ratio0 > ratio1) return 1;
        else throw std::logic_error("displacement through just on vertex");
    }
    else if(rhs.a < 0e0 && rhs.b > 0e0 && rhs.a + rhs.b > 1e0)
    {
        const T ratio1 = cross_ratio(lhs, disp, 1);
        const T ratio2 = cross_ratio(lhs, disp, 2);
             if(ratio1 < ratio2) return 1;
        else if(ratio1 > ratio2) return 2;
        else throw std::logic_error("displacement through just on vertex");
    }
    else
    {
        std::cerr << "lhs = " << lhs.a << ", " << lhs.b << std::endl;
        std::cerr << "rhs = " << rhs.a << ", " << rhs.b << std::endl;
        throw std::logic_error("invalid position");
    }
}

// return edge_id of the edge that crosses vector from lhs to rhs.
// if the vector lhs->rhs crosses some edges, return nearest one.
template<typename T>
std::size_t cross_edge(const ParametricPosition<T>& lhs,
        const ParametricPosition<T>& rhs, const T t = 0e0)
{
         if(is_in_face(lhs) && !is_in_face(rhs))
        return cross_edge_lhs_is_in_face_rhs_is_out_of_face(lhs, rhs);
    else if(!is_in_face(lhs) && is_in_face(rhs))
        return cross_edge_lhs_is_in_face_rhs_is_out_of_face(rhs, lhs);
    else
        throw std::invalid_argument("the implementation of this case is hell");
}

// XXX: NOTE: pos is relative position to the vtxs[0].
template<typename coordT>
ParametricPosition<typename value_type_helper<coordT>::type>
to_parametric(const coordT& pos, const coordT& a_vec, const coordT& b_vec,
     const typename value_type_helper<coordT>::type tol = 1e-12)
{
    typedef circular_iteration<3> triangular;
    typedef typename value_type_helper<coordT>::type valueT;

    // to solve (a * a_vec + b * b_vec == pos), compute determinant of 3 matrices
    boost::array<valueT, 3> determinants;
    determinants[0] = a_vec[0] * b_vec[1] - b_vec[0] * a_vec[1];
    determinants[1] = a_vec[1] * b_vec[2] - b_vec[1] * a_vec[2];
    determinants[2] = a_vec[2] * b_vec[0] - b_vec[2] * a_vec[0];

    // use matrix that has largest determinant
    const typename boost::array<valueT, 3>::iterator max_iter =
        std::max_element(determinants.begin(), determinants.end());
    const std::size_t max_index = std::distance(determinants.begin(), max_iter);
    assert(0 <= max_index && max_index < 3);
    if(*max_iter == 0e0) throw std::logic_error("division by zero");

    const valueT alpha =
        (b_vec[triangular::increment(max_index)] * pos[max_index] - 
         b_vec[max_index] * pos[triangular::increment(max_index)]) /
        determinants[max_index];

    const valueT beta =
        (a_vec[max_index] * pos[triangular::increment(max_index)] - 
         a_vec[triangular::increment(max_index)] * pos[max_index]) /
        determinants[max_index];

    // confirm (alpha, beta) * (a, b) == pos
    if(std::abs(alpha * a_vec[triangular::decrement(max_index)] +
                beta  * b_vec[triangular::decrement(max_index)] -
                pos[triangular::decrement(max_index)]) > tol)
    {
        throw std::logic_error("to_parametric: invalid solution");
    }

    return ParametricPosition<valueT>(alpha, beta);
}

template<typename coordT>
ParametricPosition<typename value_type_helper<coordT>::type>
to_parametric(const coordT& pos, const boost::array<coordT,3>& vtxs,
    const typename value_type_helper<coordT>::type tol = 1e-12)
{
    // make parametric axes
    const coordT a_vec = vtxs[1] - vtxs[0];
    const coordT b_vec = vtxs[2] - vtxs[0];
    return to_parametric(pos, a_vec, b_vec, tol);
}


// return value is relative position to the origin
template<typename coordT>
inline coordT
to_absolute(const ParametricPosition<typename value_type_helper<coordT>::type>& para,
            const coordT& a_vec, const coordT& b_vec)
{
    return a_vec * para.a + b_vec * para.b;
}

template<typename coordT>
inline coordT
to_absolute(const ParametricPosition<typename value_type_helper<coordT>::type>& para,
            const boost::array<coordT,3>& vtxs)
{
    const coordT a_vec  = vtxs[1] - vtxs[0];
    const coordT b_vec  = vtxs[2] - vtxs[0];
    return a_vec * para.a + b_vec * para.b;
}


template<typename coordT>
ParametricPosition<typename value_type_helper<coordT>::type>
projection(const coordT& pos,
    const boost::array<coordT,3>& vtxs,
    const coordT& normal,
    const typename value_type_helper<coordT>::type tol = 1e-12)
{
    typedef typename value_type_helper<coordT>::type valueT;

    const coordT a_vec  = vtxs[1] - vtxs[0];
    const coordT b_vec  = vtxs[2] - vtxs[0];
    const valueT det_term1 = a_vec[0] * b_vec[1] * normal[2] +
                        a_vec[1] * b_vec[2] * normal[0] +
                        a_vec[2] * b_vec[0] * normal[1];
    const valueT det_term2 = a_vec[0] * b_vec[2] * normal[1] +
                        a_vec[1] * b_vec[0] * normal[2] +
                        a_vec[2] * b_vec[1] * normal[0];
    const valueT determinant = det_term1 - det_term2;
    if(determinant == 0e0)
        throw std::logic_error("division by zero");

    const valueT alpha = ((b_vec[1]*normal[2] - b_vec[2]*normal[1])*pos[0] +
                          (b_vec[2]*normal[0] - b_vec[0]*normal[2])*pos[1] +
                          (b_vec[0]*normal[1] - b_vec[1]*normal[0])*pos[2])
                          / determinant;
    const valueT beta  = ((a_vec[2]*normal[1] - a_vec[1]*normal[2])*pos[0] +
                          (a_vec[0]*normal[2] - a_vec[2]*normal[0])*pos[1] +
                          (a_vec[1]*normal[0] - a_vec[0]*normal[1])*pos[2])
                          / determinant;
    const valueT gamma = ((a_vec[1]*b_vec[2] - a_vec[2]*b_vec[1])*pos[0] +
                          (a_vec[2]*b_vec[0] - a_vec[0]*b_vec[2])*pos[1] +
                          (a_vec[0]*b_vec[1] - a_vec[1]*b_vec[0])*pos[2])
                          / determinant;
    if(gamma > tol)
        throw std::invalid_argument("projection: over the tolerance");

    return ParametricPosition<valueT>(alpha, beta);
}

// operator+-* {{{
template<typename T>
inline ParametricPosition<T> operator+(const ParametricPosition<T>& lhs,
                                const ParametricPosition<T>& rhs)
{
    return ParametricPosition<T>(lhs.a + rhs.a, lhs.b + rhs.b);
}

template<typename T>
inline ParametricPosition<T> operator-(const ParametricPosition<T>& lhs,
                                const ParametricPosition<T>& rhs)
{
    return ParametricPosition<T>(lhs.a - rhs.a, lhs.b - rhs.b);
}

template<typename T>
inline ParametricPosition<T> operator*(const ParametricPosition<T>& lhs,
                               const T& rhs)
{
    return ParametricPosition<T>(lhs.a * rhs, lhs.b * rhs);
}

template<typename T>
inline ParametricPosition<T> operator*(const T& lhs,
                               const ParametricPosition<T>& rhs)
{
    return ParametricPosition<T>(lhs * rhs.a, lhs * rhs.b);
}
// }}}


#endif /* GFRD_POLYGON_PARAMETRIC_POSITION */
