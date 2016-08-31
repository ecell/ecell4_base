#ifndef EGFRD_POLYGON
#define EGFRD_POLYGON

#include <ecell4/core/Shape.hpp>
#include "exceptions.hpp"
#include "FaceTriangle.hpp"

struct Polygon : public ecell4::Shape
{
    Polygon(){}
    ~Polygon(){}

    // nearest
    // (intruders, (idx, (distance, max_radius)))
    std::pair<std::vector<std::size_t>, std::pair<std::size_t, std::pair<Real, Real> > >
    get_intruder_faces(const Real3& pos, const Real range) const;

    std::vector<FaceTrialgne<Real3> > faces;

// shape
    dimensiton_kind dimension() const {return THREE;}
    Real  is_inside(const Real3& coord) const;
    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const;
    bool  test_AABB(const Real3& l, const Real3& u) const;
};

#endif //EGFRD_POLYGON
