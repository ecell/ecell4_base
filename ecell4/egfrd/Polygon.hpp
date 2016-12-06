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
    typedef std::size_t face_id_type;

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

    // assume the segment collides a face at first and
    // return pairof(new begin, new end) and id of collided face
    std::pair<std::pair<coordinate_type, coordinate_type>, face_id_type>
    apply_reflection(const coordinate_type& pos, const coordinate_type& displacement,
                     const std::vector<face_id_type>& intruder_faces,
                     const face_id_type ignore_face) const;

    std::pair<std::pair<coordinate_type, coordinate_type>, face_id_type>
    apply_reflection(const coordinate_type& pos, const coordinate_type& displacement,
                     const face_id_type intruder_face) const;

    std::pair<bool, std::pair<Real, face_id_type> >
    intersect_ray(const coordinate_type& pos, const coordinate_type& disp,
                  const face_id_type ignore_face) const;

    static face_id_type make_nonsence_id(){return std::numeric_limits<std::size_t>::max();}

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
std::pair<std::pair<typename Polygon<coordT>::coordinate_type,
                    typename Polygon<coordT>::coordinate_type>,
          typename Polygon<coordT>::face_id_type>
Polygon<coordT>::apply_reflection(
        const coordinate_type& pos, const coordinate_type& displacement,
        const std::vector<face_id_type>& intruder_faces,
        const face_id_type ignore_face) const
{
    const coordinate_type end = pos + displacement;
    if(intruder_faces.empty())
        return std::make_pair(std::make_pair(end, end), make_nonsence_id());

    bool collide_face = false;
    coordinate_type next_begin = end;
    face_id_type first_collide_face_idx = make_nonsence_id();
    Real         first_collide_distance = length(displacement);

    for(typename std::vector<face_id_type>::const_iterator
        iter = intruder_faces.begin(); iter != intruder_faces.end(); ++iter)
    {
        if(*iter == ignore_face) continue;

        const std::pair<bool, coordinate_type> test_result =
            test_intersect_segment_triangle(pos, end, faces.at(*iter));

        if(test_result.first)
        {
            const Real dist_to_face = length(test_result.second - pos);

            if(dist_to_face < first_collide_distance)
            {
                collide_face = true;
                first_collide_face_idx = *iter;
                first_collide_distance = dist_to_face;
                next_begin = test_result.second;
            }
            else if(dist_to_face == first_collide_distance)
            {
                throw ecell4::NotImplemented("collide 2 object at the same time");
            }
        }
    }
    if(!collide_face)
        return std::make_pair(std::make_pair(end, end), make_nonsence_id());

    const coordinate_type next_end =
        reflect_plane(pos, end, faces.at(first_collide_face_idx));

    return std::make_pair(
            std::make_pair(next_begin, next_end), first_collide_face_idx);
}


template<typename coordT>
std::pair<std::pair<typename Polygon<coordT>::coordinate_type,
                    typename Polygon<coordT>::coordinate_type>,
          typename Polygon<coordT>::face_id_type>
Polygon<coordT>::apply_reflection(
        const coordinate_type& pos, const coordinate_type& displacement,
        const face_id_type intruder_face) const
{
    const coordinate_type end = pos + displacement;
    const std::pair<bool, coordinate_type> test_result =
        test_intersect_segment_triangle(pos, end, faces.at(intruder_face));

    const coordinate_type next_end =
        reflect_plane(pos, end, faces.at(intruder_face));

    return std::make_pair(
            std::make_pair(test_result.second, next_end), intruder_face);
}


template<typename coordT>
std::pair<bool, std::pair<Real, typename Polygon<coordT>::face_id_type> >
Polygon<coordT>::intersect_ray(
        const coordinate_type& pos, const coordinate_type& disp, 
        const face_id_type ignore_face) const
{
    const std::pair<std::vector<face_id_type>,
                    std::pair<face_id_type, std::pair<Real, Real> > > intruders =
            this->get_faces_within_radius(pos, length(disp));

    bool collide_face = false;
    face_id_type first_collide_face_idx = std::numeric_limits<std::size_t>::max();
    Real         first_collide_distance = length(disp);
    const coordinate_type end = pos + disp;
    for(typename std::vector<face_id_type>::const_iterator
        iter = intruders.first.begin(); iter != intruders.first.end(); ++iter)
    {
        if(*iter == ignore_face) continue;

        const std::pair<bool, coordinate_type> test_result =
            test_intersect_segment_triangle(pos, end, this->faces.at(*iter));

        if(test_result.first)
        {
            const Real dist_to_face = length(test_result.second - pos);
            if(dist_to_face < first_collide_distance)
            {
                collide_face = true;
                first_collide_face_idx = *iter;
                first_collide_distance = dist_to_face;
            }
        }
    }

    return std::make_pair(collide_face,
                std::make_pair(first_collide_distance, first_collide_face_idx));
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
