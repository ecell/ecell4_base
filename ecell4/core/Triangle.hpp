#ifndef ECELL_CORE_TRIANGLE
#define ECELL_CORE_TRIANGLE

#include "Shape.hpp"
#include "exceptions.hpp"

namespace ecell4
{

struct Triangle : public Shape
{
  public:

    Triangle();
    explicit Triangle(const boost::array<Real3, 3>& vertices);
    Triangle(const Real3& a, const Real3& b, const Real3& c);

    Real3 const& normal() const
    {
        return normal_;
    }
    Real3 const& represent() const
    {
        return edges_[0];
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

    boost::array<Real3, 3> const& vertices() const
    {
        return vertices_;
    }
    boost::array<Real3, 3> const& edges() const
    {
        return edges_;
    }
    boost::array<Real, 3> const& lengths_of_edges() const
    {
        return lengths_;
    }
    Real angle(const Real3& a, const Real3& b) const
    {
        return acos(dot_product(a, b) / std::sqrt(length_sq(a) * length_sq(b)));
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
    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
    {
        Real a1 = rng.uniform(0.0, 1.0);
        Real a2 = rng.uniform(0.0, 1.0);
        if(a1 + a2 > 1.0)
        {
            a1 = 1.0 - a1;
            a2 = 1.0 - a2;
        }
        return edges_[0] * a1 - edges_[2] * a2;
    }
    bool test_AABB(const Real3& l, const Real3& u) const
    {
        throw NotImplemented("Triangle::test_AABB(l, u)");
    }
    void bounding_box(
            const Real3& edge_lengths, Real3& lower, Real3& upper) const
    {
        throw NotImplemented("Triangle::bounding_box(edge_length, lower, upper)");
    }

  private:

    Real3 normal_;
    boost::array<Real, 3>  lengths_;
    boost::array<Real, 3>  angles_;
    boost::array<Real3, 3> vertices_;
    boost::array<Real3, 3> edges_;
};

}

#endif /*ECELL_CORE_TRIANGLE*/
