#ifndef ECELL4_EGFRD_FACE_TRIANGLE_HPP
#define ECELL4_EGFRD_FACE_TRIANGLE_HPP
#include <ecell4/core/geometry.hpp>
#include <array>

namespace ecell4
{
namespace egfrd
{
template<typename coordT>
struct FaceTriangle
{
  public:
    typedef coordT                                   position_type;
    typedef position_type                            vector_type;
    typedef typename element_type_of<coordT>::type length_type;
    typedef std::size_t                              size_type;
    typedef size_type                                index_type;
    typedef std::array<length_type, 3>             length_container_type;
    typedef std::array<position_type, 3>           position_container_type;

  public:
    FaceTriangle(){}
    explicit FaceTriangle(const position_container_type& vertices)
        : normal_(cross_product(vertices[1] - vertices[0],
                                vertices[2] - vertices[0]) /
                  length(cross_product(vertices[1] - vertices[0],
                                       vertices[2] - vertices[0]))),
          vertices_(vertices)
    {
        edges_[0] = vertices[1] - vertices[0];
        edges_[1] = vertices[2] - vertices[1];
        edges_[2] = vertices[0] - vertices[2];
        lengths_[0] = length(edges_[0]);
        lengths_[1] = length(edges_[1]);
        lengths_[2] = length(edges_[2]);
        angles_[0] = calc_angle(edges_[0], edges_[2] * -1.0);
        angles_[1] = calc_angle(edges_[1], edges_[0] * -1.0);
        angles_[2] = calc_angle(edges_[2], edges_[1] * -1.0);
    }

    FaceTriangle(const position_type& a, const position_type& b,
                          const position_type& c)
        : normal_(cross_product(b - a, c - a) / length(cross_product(b - a, c - a)))
    {
        vertices_[0] = a;
        vertices_[1] = b;
        vertices_[2] = c;
        edges_[0] = vertices_[1] - vertices_[0];
        edges_[1] = vertices_[2] - vertices_[1];
        edges_[2] = vertices_[0] - vertices_[2];
        lengths_[0] = length(edges_[0]);
        lengths_[1] = length(edges_[1]);
        lengths_[2] = length(edges_[2]);
        angles_[0] = calc_angle(edges_[0], edges_[2] * -1.0);
        angles_[1] = calc_angle(edges_[1], edges_[0] * -1.0);
        angles_[2] = calc_angle(edges_[2], edges_[1] * -1.0);
    }

    vector_type   const& normal()    const {return normal_;}
    vector_type   const& represent() const {return edges_[0];}
    position_type const& vertex_at        (const index_type i) const {return vertices_.at(i);}
    vector_type   const& edge_at          (const index_type i) const {return edges_.at(i);}
    length_type   const& length_of_edge_at(const index_type i) const {return lengths_.at(i);}
    length_type   const& angle_at         (const index_type i) const {return angles_.at(i);}

    position_container_type const& vertices()         const {return vertices_;}
    position_container_type const& edges()            const {return edges_;}
    length_container_type   const& lengths_of_edges() const {return lengths_;}

  private:

    vector_type             normal_;
    length_container_type   lengths_;
    length_container_type   angles_;
    position_container_type vertices_;
    position_container_type edges_;
};

namespace detail
{
template<typename coordT>
coordT closest_point(const coordT& pos, const std::array<coordT, 3>& vertices)
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
    {
        return a;
    }

    const coordT bp = pos - b;
    const valueT d3 = dot_product(ab, bp);
    const valueT d4 = dot_product(ac, bp);
    if (d3 >= 0.0 && d4 <= d3)
    {
        return b;
    }

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
    {
        return c;
    }

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
std::pair<bool, coordT>
test_intersect_segment_triangle(const coordT& begin, const coordT& end,
                                const std::array<coordT, 3>& vertices)
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
} // detail

template<typename coordT>
typename element_type_of<coordT>::type
distance_point_to_triangle(const coordT& pos, const FaceTriangle<coordT>& face)
{
    std::array<coordT, 3> triangle = face.vertices();
    if(dot_product(pos - face.vertex_at(0), face.normal()) < 0)
    {
        triangle[0] = face.vertex_at(2);
        triangle[1] = face.vertex_at(1);
        triangle[2] = face.vertex_at(0);
    }
    return length(detail::closest_point(pos, triangle) - pos);
}

template<typename coordT>
std::pair<bool, coordT>
test_intersect_segment_triangle(const coordT& begin, const coordT& end,
                                const FaceTriangle<coordT>& face)
{
    const coordT line = end - begin;
    if(dot_product(line, face.normal()) < 0.0)
    {
        return detail::test_intersect_segment_triangle(begin, end, face.vertices());
    }
    else
    {
        std::array<coordT, 3> rev;
        rev[0] = face.vertex_at(2);
        rev[1] = face.vertex_at(1);
        rev[2] = face.vertex_at(0);
        return detail::test_intersect_segment_triangle(begin, end, rev);
    }
}

template<typename coordT>
coordT reflect_plane(const coordT& begin, const coordT& end,
                     const FaceTriangle<coordT>& face)
{
    return reflect_plane(begin, end, face.normal(), face.vertex_at(0));
}

} // egfrd
} // ecell4
#endif /* GFRD_POLYGON_FACE_TRIANGLE */
