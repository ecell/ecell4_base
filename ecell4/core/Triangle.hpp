#ifndef ECELL_CORE_TRIANGLE
#define ECELL_CORE_TRIANGLE

#include "Shape.hpp"
#include "geometry.hpp"
#include "exceptions.hpp"
#include "TriangleView.hpp"
#include "BoundaryCondition.hpp"

namespace ecell4
{

struct Triangle : public Shape
{
  public:

    Triangle();
    explicit Triangle(const std::array<Real3, 3>& vertices);
    explicit Triangle(const TriangleView& tv);
    explicit Triangle(const TriangleConstView& tv);
    Triangle(const Real3& a, const Real3& b, const Real3& c);

    Triangle& operator=(const Triangle& rhs);
    Triangle& operator=(const TriangleView& tv);
    Triangle& operator=(const TriangleConstView& tv);

    Real3 const& normal() const
    {
        return normal_;
    }
    Real3 represent() const
    {
        return edges_[0] / lengths_[0];
    }

    Real3 const& vertex_at(const std::size_t i) const
    {
        return vertices_.at(i);
    }
    Real3 const& edge_at(const std::size_t i) const
    {
        return edges_.at(i);
    }
    Real const& length_of_edge_at(const std::size_t i) const
    {
        return lengths_.at(i);
    }
    Real const& angle_at(const std::size_t i) const
    {
        return angles_.at(i);
    }

    std::array<Real3, 3> const& vertices() const
    {
        return vertices_;
    }
    std::array<Real3, 3> const& edges() const
    {
        return edges_;
    }
    std::array<Real, 3> const& lengths_of_edges() const
    {
        return lengths_;
    }
    Real area() const
    {
        return 0.5 * length(cross_product(edges_[0], edges_[2] * (-1.0)));
    }

    // shape
    dimension_kind dimension() const
    {
        return TWO;
    }
    Real is_inside(const Real3& coord) const
    {
        throw NotImplemented("Triangle::is_inside(coord)");
    }
    Real3 draw_position(std::shared_ptr<RandomNumberGenerator>& rng) const
    {
        Real a1 = rng->uniform(0.0, 1.0);
        Real a2 = rng->uniform(0.0, 1.0);
        if(a1 + a2 > 1.0)
        {
            a1 = 1.0 - a1;
            a2 = 1.0 - a2;
        }
        return vertices_[0] + edges_[0] * a1 - edges_[2] * a2;
    }
    bool test_AABB(const Real3& l, const Real3& u) const
    {
        throw NotImplemented("Triangle::test_AABB(l, u)");
    }
    void bounding_box(
            const Real3& edge_lengths, Real3& lower, Real3& upper) const
    {
        const auto& v1 = this->vertices_[0];
        const auto& v2 = this->vertices_[1];
        const auto& v3 = this->vertices_[2];

        lower[0] = std::max(0.0, std::min(std::min(v1[0], v2[0]), v3[0]));
        lower[1] = std::max(0.0, std::min(std::min(v1[1], v2[1]), v3[1]));
        lower[2] = std::max(0.0, std::min(std::min(v1[2], v2[2]), v3[2]));

        upper[0] = std::min(edge_lengths[0], std::max(std::max(v1[0], v2[0]), v3[0]));
        upper[1] = std::min(edge_lengths[1], std::max(std::max(v1[1], v2[1]), v3[1]));
        upper[2] = std::min(edge_lengths[2], std::max(std::max(v1[2], v2[2]), v3[2]));
        return;
    }

  private:

    Real3 normal_;
    std::array<Real, 3>  lengths_;
    std::array<Real, 3>  angles_;
    std::array<Real3, 3> vertices_;
    std::array<Real3, 3> edges_;
};

Real distance_sq_point_Triangle(const Real3& pos, const Triangle& tri);

inline Real distance_point_Triangle(const Real3& pos, const Triangle& tri)
{
    return std::sqrt(distance_sq_point_Triangle(pos, tri));
}

// ------------------------------------------------------------------------
// considering PBC

Real distance_sq_point_Triangle(const Real3& pos, const Triangle& tri,
                                const Boundary& b);
Real distance_sq_point_Triangle(const Real3& pos, const Triangle& tri,
                                const std::unique_ptr<Boundary>& b);
Real distance_sq_point_Triangle(const Real3& pos, const Triangle& tri,
                                const std::shared_ptr<Boundary>& b);

inline Real distance_point_Triangle(const Real3& pos, const Triangle& tri,
                                    const Boundary& b)
{
    return std::sqrt(distance_sq_point_Triangle(pos, tri, b));
}
inline Real distance_point_Triangle(const Real3& pos, const Triangle& tri,
                                    const std::unique_ptr<Boundary>& b)
{
    return std::sqrt(distance_sq_point_Triangle(pos, tri, b));
}
inline Real distance_point_Triangle(const Real3& pos, const Triangle& tri,
                                    const std::shared_ptr<Boundary>& b)
{
    return std::sqrt(distance_sq_point_Triangle(pos, tri, b));
}

template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Triangle& tri)
{
    os << "Triangle(" << tri.vertex_at(0) << ", " << tri.vertex_at(1) << ", "
       << tri.vertex_at(2) << ')';
    return os;
}

} // ecell
#endif /*ECELL_CORE_TRIANGLE*/
