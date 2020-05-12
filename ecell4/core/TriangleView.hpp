#ifndef ECELL_CORE_TRIANGLE_VIEW
#define ECELL_CORE_TRIANGLE_VIEW

#include "Shape.hpp"
#include "geometry.hpp"
#include "exceptions.hpp"

namespace ecell4
{

struct TriangleView : public Shape
{
  public:

    TriangleView(Real3& a, Real3& b, Real3& c)
    {
        this->vertices_[0] = &a;
        this->vertices_[1] = &b;
        this->vertices_[2] = &c;
    }

    TriangleView(TriangleView const& rhs)
    {
        this->vertices_[0] = rhs.vertices_[0];
        this->vertices_[1] = rhs.vertices_[1];
        this->vertices_[2] = rhs.vertices_[2];
    }
    TriangleView& operator=(TriangleView const& rhs)
    {
        this->vertices_[0] = rhs.vertices_[0];
        this->vertices_[1] = rhs.vertices_[1];
        this->vertices_[2] = rhs.vertices_[2];
        return *this;
    }

    Real3 normal() const
    {
        const Real3 n = cross_product(this->edges(0), this->edges(2) * (-1));
        return (n * (1. / length(n)));
    }
    Real3 represent() const
    {
        const Real3 e = this->edges(0);
        return e / length(e);
    }

    Real3& vertex_at(const std::size_t i) const
    {
        return *(vertices_.at(i));
    }
    Real3& vertices(const std::size_t i) const
    {
        return *(vertices_[i]);
    }

    Real3 edge_at(const std::size_t i) const
    {
        return *(vertices_.at(i==2?0:i+1)) - *(vertices_.at(i));
    }
    Real3 edges(const std::size_t i) const
    {
        return *(vertices_[i==2?0:i+1]) - *(vertices_[i]);
    }

    Real length_of_edge_at(const std::size_t i) const
    {
        return length(edge_at(i));
    }

    Real angle_at(const std::size_t i) const
    {
        return calc_angle(this->edges(i), this->edges(i==0?2:i-1) * -1.0);
    }

    Real area() const
    {
        return 0.5 * length(cross_product(this->edges(0), this->edges(2) * (-1.0)));
    }

    // shape
    dimension_kind dimension() const
    {
        return TWO;
    }
    Real is_inside(const Real3& coord) const
    {
        throw NotImplemented("TriangleView::is_inside(coord)");
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
        return *(this->vertices_[0]) + this->edges(0) * a1 - this->edges(2) * a2;
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

    std::array<Real3*, 3> const& get_vertices() const {return vertices_;}

  private:
    std::array<Real3*, 3> vertices_;
};

struct TriangleConstView : public Shape
{
  public:

    TriangleConstView(const Real3& a, const Real3& b, const Real3& c)
    {
        this->vertices_[0] = &a;
        this->vertices_[1] = &b;
        this->vertices_[2] = &c;
    }

    TriangleConstView(TriangleConstView const& rhs)
    {
        this->vertices_[0] = rhs.vertices_[0];
        this->vertices_[1] = rhs.vertices_[1];
        this->vertices_[2] = rhs.vertices_[2];
    }
    TriangleConstView& operator=(TriangleConstView const& rhs)
    {
        this->vertices_[0] = rhs.vertices_[0];
        this->vertices_[1] = rhs.vertices_[1];
        this->vertices_[2] = rhs.vertices_[2];
        return *this;
    }
    TriangleConstView(TriangleView const& rhs)
    {
        this->vertices_[0] = rhs.get_vertices()[0];
        this->vertices_[1] = rhs.get_vertices()[1];
        this->vertices_[2] = rhs.get_vertices()[2];
    }
    TriangleConstView& operator=(TriangleView const& rhs)
    {
        this->vertices_[0] = rhs.get_vertices()[0];
        this->vertices_[1] = rhs.get_vertices()[1];
        this->vertices_[2] = rhs.get_vertices()[2];
        return *this;
    }

    Real3 normal() const
    {
        const Real3 n = cross_product(this->edges(0), this->edges(2) * (-1));
        return (n * (1. / length(n)));
    }
    Real3 represent() const
    {
        const Real3 e = this->edges(0);
        return e / length(e);
    }

    Real3 const& vertex_at(const std::size_t i) const
    {
        return *(vertices_.at(i));
    }
    Real3 const& vertices(const std::size_t i) const
    {
        return *(vertices_[i]);
    }

    Real3 edge_at(const std::size_t i) const
    {
        return *(vertices_.at(i==2?0:i+1)) - *(vertices_.at(i));
    }
    Real3 edges(const std::size_t i) const
    {
        return *(vertices_[i==2?0:i+1]) - *(vertices_[i]);
    }

    Real length_of_edge_at(const std::size_t i) const
    {
        return length(edge_at(i));
    }

    Real angle_at(const std::size_t i) const
    {
        return calc_angle(this->edges(i), this->edges(i==0?2:i-1) * -1.0);
    }

    Real area() const
    {
        return 0.5 * length(cross_product(this->edges(0), this->edges(2) * (-1.0)));
    }

    // shape
    dimension_kind dimension() const
    {
        return TWO;
    }
    Real is_inside(const Real3& coord) const
    {
        throw NotImplemented("TriangleView::is_inside(coord)");
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
        return *(this->vertices_[0]) + this->edges(0) * a1 - this->edges(2) * a2;
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

    std::array<const Real3*, 3> const& get_vertices() const {return vertices_;}

  private:
    std::array<const Real3*, 3> vertices_;
};


} // ecell4

#endif /*ECELL_CORE_TRIANGLE_VIEW*/
