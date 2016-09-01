#include "Polygon.hpp"
#include "TriangleOperation.hpp"
#include <limits>

std::pair<std::vector<std::size_t>, std::pair<std::size_t, std::pair<Real, Real> > >
Polygon::get_faces_within_radius(const coordinate_type& pos, const Real range) const
{
    std::vector<std::size_t> intruders;
    std::size_t nearest_idx = std::numeric_limits<std::size_t>::max();
    std::pair<Real, Real> min_dist(std::numeric_limits<Real>::max(), 0.);
    std::size_t idx = 0;
    for(std::vector<FaceTriangle<coordinate_type> >::const_iterator
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

Polygon::coordinate_type
Polygon::apply_reflection(const coordinate_type& pos, const coordinate_type& displacement,
                          const std::vector<std::size_t>& intruder_idxs) const
{
    if(intruder_idxs.empty())
    {
        return pos + displacement;
    }
    else
    {
        coordinate_type begin = pos;
        coordinate_type end = pos + displacement;
        std::vector<std::pair<std::size_t, Real> > pierce;
        pierce.reserve(intruder_idxs.size());
        while(true)
        {
            pierce.clear();
            for(std::vector<std::size_t>::const_iterator
                iter = intruder_idxs.begin(); iter != intruder_idxs.end(); ++iter)
            {
                const std::pair<bool, coordinate_type> dist =
                    is_pierce(begin, end, faces.at(*iter));
                if(dist.first)
                    pierce.push_back(std::make_pair(*iter, length(dist.second)));
            }
            if(pierce.empty()) break; // particle goes through no faces. return.

            // check which face the particle collides first.
            std::vector<std::size_t> first_collision_index_list;
            Real minimum_distance = std::numeric_limits<Real>::max();
            std::size_t multiplicity = 0;
            for(std::vector<std::pair<std::size_t, Real> >::const_iterator
                iter = pierce.begin(); iter != pierce.end(); ++iter)
            {
                if(iter->second == minimum_distance)
                {
                    first_collision_index_list.push_back(iter->first);
                }
                else if(iter->second < minimum_distance)
                {
                    first_collision_index_list.clear();
                    first_collision_index_list.push_back(iter->first);
                }
            }

            // update begin and end
            if(first_collision_index_list.size() == 1)
            {
                coordinate_type tmp_begin = is_pierce(begin, end,
                        faces.at(first_collision_index_list.first())).second;
                end = reflect_plane(begin, end,
                        faces.at(first_collision_index_list.first()));
                begin = tmp_begin;
            }
            else if(first_collision_index_list.size() > 1)
            {
                throw ecell4::NotImplemented("particle collides several faces at once");
            }
            else
            {
                throw std::logic_error("app_reflect: never reach here");
            }
        }
        return end;
    }
}

Real Polygon::is_inside(const coordinate_type& coord) const
{
    throw ecell4::NotImplemented("polygon::is_inside");
}

Polygon::coordinate_type
Polygon::draw_position(boost::shared_ptr<ecell4::RandomNumberGenerator>& rng) const
{
    throw ecell4::NotImplemented("polygon::draw_position");
}

bool Polygon::test_AABB(const coordinate_type& l, const coordinate_type& u) const
{
    throw ecell4::NotImplemented("polygon::test_AABB");
}
