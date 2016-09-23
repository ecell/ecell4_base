#ifndef EGFRD_POLYGON
#define EGFRD_POLYGON

#include <ecell4/core/Shape.hpp>
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
        const std::vector<std::size_t>& intruder_idxs) const
{
    if(intruder_idxs.empty())
    {
        return pos + displacement;
    }
    else
    {
//         std::cerr << "length of disp " << length(displacement) << std::endl;
//         std::cerr << "size of intruders  " << intruder_idxs.size() << std::endl;

        coordinate_type begin = pos;
        coordinate_type end = pos + displacement;
        std::vector<std::pair<std::size_t, Real> > pierce;
        pierce.reserve(intruder_idxs.size());
        std::size_t loop_time = 0;
        std::size_t ignore_idx = std::numeric_limits<std::size_t>::max();
        while(true)
        {
            ++loop_time;

            pierce.clear();
            for(std::vector<std::size_t>::const_iterator
                iter = intruder_idxs.begin(); iter != intruder_idxs.end(); ++iter)
            {
                if(*iter == ignore_idx) continue;
                const std::pair<bool, coordinate_type> dist =
                    is_pierce(begin, end, faces.at(*iter));
                if(dist.first)
                    pierce.push_back(std::make_pair(*iter, length(dist.second - begin)));
            }
            if(pierce.empty()) break; // particle goes through no faces. return.

            // check which face the particle collides first.
            std::vector<std::size_t> first_collision_index_list;
            Real minimum_distance = std::numeric_limits<Real>::max();
            for(std::vector<std::pair<std::size_t, Real> >::const_iterator
                iter = pierce.begin(); iter != pierce.end(); ++iter)
            {
                if(iter->second == minimum_distance)
                {
                    first_collision_index_list.push_back(iter->first);
                }
                else if(iter->second < minimum_distance)
                {
                    minimum_distance = iter->second;
                    first_collision_index_list.clear();
                    first_collision_index_list.push_back(iter->first);
                }
            }

//             {
//                 std::cerr << "now loop " << loop_time << std::endl;
//                 std::cerr << "begin " << begin << std::endl;
//                 std::cerr << "end   " << end   << std::endl;
//                 std::cerr << "pierce.size " << pierce.size() << std::endl;
//             }
            // update begin and end
            if(first_collision_index_list.size() == 1)
            {
                const coordinate_type tmp_begin = is_pierce(begin, end,
                        faces.at(first_collision_index_list.front())).second;
                const coordinate_type tmp_end = reflect_plane(begin, end,
                        faces.at(first_collision_index_list.front()));
                assert(length(tmp_end - end) > 1e-12);
//                 if(length(tmp_end - end) < 1e-12)
//                 {
//                     const std::pair<bool, coordinate_type> d =
//                         is_pierce(begin, end, faces.at(first_collision_index_list.front()));
//                     std::cerr << tmp_begin << std::endl;
//                     std::cerr << tmp_end << std::endl;
//                     std::cerr << d.first << ", " << d.second << std::endl;
//                     std::cerr << length(tmp_end - end) << std::endl;
//                     assert(false);
//                 }
                end = tmp_end;
                begin = tmp_begin;
                ignore_idx = first_collision_index_list.front();
            }
            else if(first_collision_index_list.size() > 1)
            {
                throw ecell4::NotImplemented("particle collides several faces at once");
            }
            else
            {
                throw std::logic_error("app_reflect: never reach here");
            }

//             if(loop_time > 10)
//                 throw std::logic_error("too much reflection");
        }
        return end;
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
