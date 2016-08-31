#include "Polygon.hpp"
#include "TriangleOperation.hpp"
#include <limits>

std::pair<std::vector<std::size_t>, std::pair<std::size_t, std::pair<Real, Real> > >
Polygon::get_intruder_faces(const Real3& pos, const Real range) const
{
    std::vector<std::size_t> intruders;
    std::size_t nearest_idx = std::numeric_limits<std::size_t>::max();
    std::pair<Real, Real> min_dist(std::numeric_limits<Real>::max(), 0.);
    std::size_t idx = 0;
    for(typename std::vector<FaceTrialgne<Real3> >::const_iterator
            iter = faces.begin(); iter != faces.end(); ++iter)
    {
        std::pair<Real, Real> dist = distance(pos, *iter);
        if(dist.first <= range) // is intruder face
            intruders.push_back(idx);

        if(dist.first < min_dist.first) // is nearest one
        {
            min_dist = dist;
            nearest_idx = idx;
        }
        ++idx;
    }
    return std::make_pair(intruders, std::make_pair(nearest_idx, min_dist));
}

Real Polygon::is_inside(const Real3& coord) const
{
    throw ecell::NotImplemented("polygon::is_inside");
}

Real3 Polygon::draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const
{
    throw ecell::NotImplemented("polygon::draw_position");
}

bool Polygon::test_AABB(const Real3& l, const Real3& u) const
{
    throw ecell::NotImplemented("polygon::test_AABB");
}
