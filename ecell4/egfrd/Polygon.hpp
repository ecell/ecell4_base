#ifndef EGFRD_POLYGON
#define EGFRD_POLYGON

#include <ecell4/core/Shape.hpp>
#include <ecell4/core/Real3.hpp>
#include "exceptions.hpp"
#include "FaceTriangle.hpp"

struct Polygon : public ecell4::Shape
{
    typedef ecell4::Real3 coordinate_type;
    typedef FaceTriangle<coordinate_type> face_type;

    Polygon(){}
    ~Polygon(){}

    // nearest
    // (intruders, (idx, (distance, max_radius)))
    std::pair<std::vector<std::size_t>, std::pair<std::size_t, std::pair<Real, Real> > >
    get_faces_within_radius(const coordinate_type& pos, const Real range) const;

    void emplace(const boost::array<coordinate_type, 3>& vertices)
    {
        this->faces.push_back(face_type(vertices));
    }

    coordinate_type
    apply_reflection(const coordinate_type& pos,
                     const coordinate_type& displacement,
                     const std::vector<std::size_t>& intruder_idxs) const;

// data member
    std::vector<face_type> faces;

// for shapes (not implemented yet)
    dimension_kind dimension() const {return THREE;}
    Real  is_inside(const coordinate_type& coord) const;
    coordinate_type draw_position(boost::shared_ptr<ecell4::RandomNumberGenerator>& rng) const;
    bool  test_AABB(const coordinate_type& l, const coordinate_type& u) const;
};

#endif //EGFRD_POLYGON
