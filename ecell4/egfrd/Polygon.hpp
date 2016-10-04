#ifndef EGFRD_POLYGON
#define EGFRD_POLYGON

#include <ecell4/core/Shape.hpp>
#include <ecell4/core/AABB.hpp>
#include <ecell4/core/Real3.hpp>
#include "exceptions.hpp"
#include "FaceTriangle.hpp"
#include <limits>
#include <algorithm>

template<typename coordT>
struct Polygon : public ecell4::Shape
{
    typedef coordT coordinate_type;
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
    apply_reflection(const coordinate_type& pos, const coordinate_type& displacement,
                     const std::vector<std::size_t>& intruder_idxs,
                     const coordinate_type& edge_lengths) const;

// data member
    std::vector<face_type> faces;

// for shapes (not implemented yet)
    dimension_kind dimension() const {return THREE;}
    Real  is_inside(const coordinate_type& coord) const;
    coordinate_type draw_position(boost::shared_ptr<ecell4::RandomNumberGenerator>& rng) const;
    bool  test_AABB(const coordinate_type& l, const coordinate_type& u) const;
};

template<typename coordT>
std::pair<std::vector<std::size_t>, std::pair<std::size_t, std::pair<Real, Real> > >
Polygon<coordT>::get_faces_within_radius(const coordinate_type& pos, const Real range) const
{
    std::vector<std::size_t> intruders;
    std::size_t nearest_idx = std::numeric_limits<std::size_t>::max();
    std::pair<Real, Real> min_dist(std::numeric_limits<Real>::max(), 0.);
    std::size_t idx = 0;
    for(typename std::vector<FaceTriangle<coordinate_type> >::const_iterator
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

template<typename coordT>
typename Polygon<coordT>::coordinate_type
Polygon<coordT>::apply_reflection(
        const coordinate_type& pos, const coordinate_type& displacement,
        const std::vector<std::size_t>& intruder_idxs,
        const coordinate_type& edge_lengths/* = unit cell of periodic boundary*/) const
{
    const ecell4::AABB unitcell(coordinate_type(0.,0.,0.), edge_lengths);
    coordinate_type begin = pos;
    coordinate_type end   = pos + displacement;
    std::pair<bool, Real> collide_world = unitcell.intersect_ray(pos, displacement);
    if(collide_world.first)
    {// second is just a parameter
        collide_world.second = collide_world.second * length(displacement);
    }
    else
    {
        collide_world.second = std::numeric_limits<Real>::infinity();
    }

    if(intruder_idxs.empty() && not collide_world.first)
        return end;

    std::vector<std::size_t> intruder_faces = intruder_idxs;
    std::size_t ignore_idx = std::numeric_limits<std::size_t>::max();

    while(true)
    {
        bool        collide_face = false;
        std::size_t first_collide_face_idx = std::numeric_limits<std::size_t>::max();
        Real        first_collide_distance = collide_world.second;// inf or dist to unit cell

        for(typename std::vector<std::size_t>::const_iterator
            iter = intruder_faces.begin(); iter != intruder_faces.end(); ++iter)
        {
            if(*iter == ignore_idx) continue;

            const std::pair<bool, coordinate_type> dist =
                is_pierce(begin, end, faces.at(*iter));
            if(dist.first)
            {
                const Real dist_to_face = length(dist.second - begin);
                if(dist_to_face == collide_world.second)
                {
                    throw ecell4::NotImplemented("collide 2 object at the same time");
                }
                else if(dist_to_face < first_collide_distance)
                {
                    collide_face = true;
                    first_collide_face_idx = *iter;
                    first_collide_distance = dist_to_face;
                }
            }
        }
        if(not collide_face && not collide_world.first) return end;

        if(not collide_face)
        {// collide unit cell first

            const coordinate_type upto_collision =
                (end - begin) * (collide_world.second / length(end - begin));
            coordinate_type move_into_unit_cell;
            coordinate_type tmp_begin = begin + upto_collision;

            // XXX omg!!!
            if(std::abs(tmp_begin[0]) < 1e-12)
                move_into_unit_cell = coordinate_type(edge_lengths[0], 0., 0.);
            else if(std::abs(tmp_begin[0] - edge_lengths[0]) < 1e-12)
                move_into_unit_cell = coordinate_type(-1. * edge_lengths[0], 0., 0.);
            else if(std::abs(tmp_begin[1]) < 1e-12)
                move_into_unit_cell = coordinate_type(0., edge_lengths[1], 0.);
            else if(std::abs(tmp_begin[1] - edge_lengths[1]) < 1e-12)
                move_into_unit_cell = coordinate_type(0., -1. * edge_lengths[1], 0.);
            else if(std::abs(tmp_begin[2]) < 1e-12)
                move_into_unit_cell = coordinate_type(0., 0., edge_lengths[2]);
            else if(std::abs(tmp_begin[2] - edge_lengths[2]) < 1e-12)
                move_into_unit_cell = coordinate_type(0., 0., -1. * edge_lengths[2]);
            else throw std::logic_error("never reach here");

            begin = (tmp_begin + move_into_unit_cell);
            end   = (end       + move_into_unit_cell);

            const coordinate_type new_displ = end - begin;
            const Real length_new_displ = length(new_displ);
            collide_world = unitcell.intersect_ray(begin, new_displ);
            if(collide_world.first)
                collide_world.second = collide_world.second * length_new_displ;
            else
                collide_world.second = std::numeric_limits<Real>::infinity();

            intruder_faces = (this->get_faces_within_radius(begin, length_new_displ)).first;
            ignore_idx = std::numeric_limits<std::size_t>::max();
        }
        else
        {// collide certain face first

            const coordinate_type tmp_begin = is_pierce(begin, end,
                    faces.at(first_collide_face_idx)).second;
            const coordinate_type tmp_end = reflect_plane(begin, end,
                    faces.at(first_collide_face_idx));

            begin = tmp_begin;
            end   = tmp_end;
            ignore_idx = first_collide_face_idx;

            const coordinate_type new_displ = end - begin;
            collide_world = unitcell.intersect_ray(begin, new_displ);
            if(collide_world.first)
                collide_world.second = collide_world.second * length(new_displ);
            else
                collide_world.second = std::numeric_limits<Real>::infinity();
        }
    }
}

template<typename coordT>
Real Polygon<coordT>::is_inside(const coordinate_type& coord) const
{
    throw ecell4::NotImplemented("polygon::is_inside");
}

template<typename coordT>
typename Polygon<coordT>::coordinate_type
Polygon<coordT>::draw_position(boost::shared_ptr<ecell4::RandomNumberGenerator>& rng) const
{
    throw ecell4::NotImplemented("polygon::draw_position");
}

template<typename coordT>
bool Polygon<coordT>::test_AABB(const coordinate_type& l, const coordinate_type& u) const
{
    throw ecell4::NotImplemented("polygon::test_AABB");
}

#endif //EGFRD_POLYGON
