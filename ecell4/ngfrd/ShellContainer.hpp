#ifndef ECELL4_NGFRD_SHELL_CONTAINER_HPP
#define ECELL4_NGFRD_SHELL_CONTAINER_HPP
#include <ecell4/core/PeriodicRTree.hpp>
#include <ecell4/core/exceptions.hpp>
#include <ecell4/ngfrd/ShellID.hpp>
#include <ecell4/ngfrd/Shell.hpp>

namespace ecell4
{
namespace ngfrd
{

struct ShellAABBGetterVisitor
    : boost::static_visitor<AABB>
{
    Real margin;
    explicit ShellAABBGetterVisitor(const Real m): margin(m) {}

    AABB operator()(const SphericalShell& sh) const noexcept
    {
        const Real3 radius(sh.shape().radius() * (1.0 + margin),
                           sh.shape().radius() * (1.0 + margin),
                           sh.shape().radius() * (1.0 + margin));
        return AABB(sh.position() - radius, sh.position() + radius);
    }
    AABB operator()(const CylindricalShell& sh) const noexcept
    {
        // Here we overestimate the size of the cylinder. but overestimation is
        // okay (of course, from the viewpoint of efficiency, its not okay).
        const auto& cyl_r = sh.shape().radius();
        const auto& cyl_h = sh.shape().half_height();
        const Real  r = std::sqrt(cyl_r * cyl_r + cyl_h * cyl_h) * (1.0 + margin);

        return AABB(sh.position() - Real3(r, r, r), sh.position() + Real3(r, r, r));
    }
    AABB operator()(const CircularShell& sh) const noexcept
    {
        // Consider the case where a shell lies on two faces. In such a case,
        // the shell is no longer planner.
        const Real  r = sh.shape().radius() * (1.0 + margin);
        return AABB(sh.position() - Real3(r, r, r), sh.position() + Real3(r, r, r));
    }
    AABB operator()(const ConicalShell& sh) const noexcept
    {
        // It looks like a stupid overestimation, but it is simpler than
        // calculating the direction of a cone and including all the edges on
        // a polyogn face that involves to the shell.
        const auto& apx = sh.shape().apex();
        const auto& slt = sh.shape().slant_height();
        const Real  r = slt * (1.0 + margin);
        return AABB(apx - Real3(r, r, r), apx + Real3(r, r, r));
    }
};

class ShellContainer
{
public:
    struct ShellAABBGetter
    {
        AABB operator()(const Shell& sh, const Real margin) const noexcept
        {
            return visit(ShellAABBGetterVisitor(margin), sh);
        }
    };

    // 3D part
    using rtree_type = PeriodicRTree<ShellID, Shell, ShellAABBGetter>;
    using box_type               = typename rtree_type::box_type;
    using value_type             = typename rtree_type::value_type;
    using key_to_value_map_type  = typename rtree_type::key_to_value_map_type;
    using container_type         = typename rtree_type::container_type;
    using iterator               = typename rtree_type::iterator;
    using const_iterator         = typename rtree_type::const_iterator;

    // 2D part (contained by another container)
    using polygon_container_type = PolygonContainer<ShellID>;

public:

    ShellContainer(const Real3& edge_lengths, std::shared_ptr<Polygon> poly,
                   const Real margin = 0.1)
        : rtree_(edge_lengths, margin), polygon_(std::move(poly))
    {}

    std::pair<ShellID, Shell> get_shell(const ShellID& sid) const
    {
        if(!rtree_.has(sid))
        {
            throw_exception<NotFound>("ngfrd::ShellContainer::get_shell: "
                    "No such shell (", sid, ").");
        }
        return rtree_.get(sid);
    }

    bool update_shell(const ShellID& sid, const Shell& sh)
    {
        const auto retval = rtree_.update(sid, sh);
        assert(rtree_.diagnosis()); // just want to make it sure
        if(sh.is_circular())
        {
            poly_con_.update(sid, sh.as_circular().fid());
        }
        else if(sh.is_conical())
        {
            poly_con_.update(sid, sh.as_conical().vid());
        }
        return retval;
    }

    bool has_shell(const ShellID& sid) const
    {
        return rtree_.has(sid);
    }

    void remove_shell(const ShellID& sid)
    {
        if(!rtree_.has(sid))
        {
            throw_exception<NotFound>("ngfrd::ShellContainer::remove_shell: "
                    "No such shell (", sid, ").");
        }
        const auto& sh = rtree_.get(sid).second;
        rtree_.erase(sid, sh);

        if(sh.is_circular())
        {
            poly_con_.remove(sid, sh.as_circular().fid());
        }
        else if(sh.is_conical())
        {
            poly_con_.remove(sid, sh.as_conical().vid());
        }
        return ;
    }

    // ------------------------------------------------------------------------
    // 2D/3D checking

    boost::optional<FaceID> on_which_face(const ShellID& sid) const noexcept
    {
        return poly_con_.on_which_face(sid);
    }
    boost::optional<VertexID> on_which_vertex(const ShellID& sid) const noexcept
    {
        return poly_con_.on_which_vertex(sid);
    }

    boost::optional<std::vector<ShellID> const&>
    shells_on(const FaceID& fid) const noexcept
    {
        return poly_con_.objects_on(fid);
    }
    boost::optional<std::vector<ShellID> const&>
    shells_on(const VertexID& vid) const noexcept
    {
        return poly_con_.objects_on(vid);
    }

    // ------------------------------------------------------------------------
    // list_shells 3D

    std::vector<std::pair<std::pair<ShellID, Shell>, Real>>
    list_shells_within_radius_3D(const Real3& pos, const Real& radius) const
    {
        return list_shells_within_radius_3D_impl(pos, radius,
            [this](const ShellID& sid) noexcept -> bool {
                return this->poly_con_.on_face(sid) ||
                       this->poly_con_.on_vertex(sid) ;
            });
    }
    std::vector<std::pair<std::pair<ShellID, Shell>, Real>>
    list_shells_within_radius_3D(const Real3& pos, const Real& radius,
            const ShellID& ignore) const
    {
        return list_shells_within_radius_3D_impl(pos, radius,
            [this, ignore](const ShellID& sid) noexcept -> bool {
                return sid == ignore                ||
                       this->poly_con_.on_face(sid) ||
                       this->poly_con_.on_vertex(sid) ;
            });
    }
    std::vector<std::pair<std::pair<ShellID, Shell>, Real>>
    list_shells_within_radius_3D(const Real3& pos, const Real& radius,
            const ShellID& ignore1, const ShellID& ignore2) const
    {
        return list_shells_within_radius_3D_impl(pos, radius,
            [this, ignore1, ignore2](const ShellID& sid) noexcept -> bool {
                return sid == ignore1 || sid == ignore2 ||
                       this->poly_con_.on_face(sid)     ||
                       this->poly_con_.on_vertex(sid) ;
            });
    }

    // ------------------------------------------------------------------------
    // 2D

    std::vector<std::pair<std::pair<ShellID, Shell>, Real>>
    list_shells_within_radius_2D(const std::pair<Real3, FaceID>& pos,
                                 const Real& radius) const
    {
        return list_shells_within_radius_2D_impl(pos, radius,
            [](const ShellID&) {return false;});
    }
    std::vector<std::pair<std::pair<ShellID, Shell>, Real>>
    list_shells_within_radius_2D(const std::pair<Real3, FaceID>& pos,
            const Real& radius, const ShellID& ignore) const
    {
        return list_shells_within_radius_2D_impl(pos, radius,
            [&](const ShellID& sid) {return sid == ignore;});
    }
    std::vector<std::pair<std::pair<ShellID, Shell>, Real>>
    list_shells_within_radius_2D(const std::pair<Real3, FaceID>& pos,
            const Real& radius, const ShellID& ignore1, const ShellID& ignore2) const
    {
        return list_shells_within_radius_2D_impl(pos, radius,
            [&](const ShellID& sid) {return sid == ignore1 || sid == ignore2;});
    }

    // ------------------------------------------------------------------------
    // container like access

    container_type const& list_shells() const noexcept {return rtree_.list_objects();}

    // do not modify the shape of a shell via this reference.
    value_type&       at(const ShellID& sid)       {return this->rtree_.at(sid);}
    value_type const& at(const ShellID& sid) const {return this->rtree_.at(sid);}

    value_type&       front()       noexcept {return this->rtree_.front();}
    value_type const& front() const noexcept {return this->rtree_.front();}

    value_type&       back()       noexcept {return this->rtree_.back();}
    value_type const& back() const noexcept {return this->rtree_.back();}

    iterator        begin()       noexcept {return this->rtree_.begin();}
    iterator        end()         noexcept {return this->rtree_.end();}
    const_iterator  begin() const noexcept {return this->rtree_.begin();}
    const_iterator  end()   const noexcept {return this->rtree_.end();}
    const_iterator cbegin() const noexcept {return this->rtree_.cbegin();}
    const_iterator cend()   const noexcept {return this->rtree_.cend();}

    // ------------------------------------------------------------------------
    // reset

    void clear()
    {
        poly_con_.clear();
        rtree_.clear();
        return;
    }

    void reset_boundary(const Real3& edges)
    {
        this->rtree_.reset_boundary(edges);
        return;
    }

    // ------------------------------------------------------------------------
    // testing

    bool diagnosis() const {return this->rtree_.diagnosis();}

private:

    // ------------------------------------------------------------------------
    // query objects

    template<typename Filter>
    struct IntersectionQuery
    {
        Real3  center;
        Real   radius;
        Filter ignores;

        IntersectionQuery(const Real3& c, const Real r, Filter f) noexcept
            : center(c), radius(r), ignores(std::move(f))
        {}

        // If it does not matches, return boost::none.
        // If it matches, return pairof(shell, distance).
        boost::optional<std::pair<value_type, Real>>
        operator()(const value_type& sidp, const PeriodicBoundary& pbc) const noexcept
        {
            if(ignores(sidp.first)){return boost::none;}

            // use the same algorithm as the ParticleSpaceVectorImpl.
            const auto rhs = pbc.periodic_transpose(sidp.second.position(),
                                                    this->center);

            const auto dist = visit(ShellDistanceCalculator(this->center, pbc),
                                    sidp.second);
            if(dist <= this->radius)
            {
                return std::make_pair(sidp, dist);
            }
            return boost::none;
        }

        bool operator()(const AABB& box, const PeriodicBoundary& pbc) const noexcept
        {
            return this->distance_sq(box, this->center, pbc) <=
                   this->radius * this->radius;
        }

        // -------------------------------------------------------------------
        // AABB-point distance under the PBC
        Real distance_sq(const AABB& box, Real3 pos, const PeriodicBoundary& pbc) const noexcept
        {
            pos = pbc.periodic_transpose(pos, (box.upper() + box.lower()) * 0.5);

            Real dist_sq = 0;
            for(std::size_t i=0; i<3; ++i)
            {
                const auto v = pos[i];
                if(v < box.lower()[i])
                {
                    dist_sq += (v - box.lower()[i]) * (v - box.lower()[i]);
                }
                else if(box.upper()[i] < v)
                {
                    dist_sq += (v - box.upper()[i]) * (v - box.upper()[i]);
                }
            }
            return dist_sq;
        }
    };

    template<typename Filter>
    std::vector<std::pair<std::pair<ShellID, Shell>, Real>>
    list_shells_within_radius_3D_impl(
            const Real3& pos, const Real& radius, Filter filter) const
    {
        std::vector<std::pair<std::pair<ShellID, Shell>, Real>> retval;

        rtree_.query(IntersectionQuery<Filter>(pos, radius, std::move(filter)),
                     std::back_inserter(retval));

        std::sort(retval.begin(), retval.end(), utils::pair_second_element_comparator<
                  std::pair<ShellID, Shell>, Real>());
        return retval;
    }

    template<typename Filter>
    std::vector<std::pair<std::pair<ShellID, Shell>, Real>>
    list_shells_within_radius_2D_impl(const std::pair<Real3, FaceID>& center,
            const Real& radius, Filter filter) const
    {
        std::vector<std::pair<std::pair<ShellID, Shell>, Real>> retval;

        const auto& pos = center.first;
        const auto& fid = center.second;

        // check shells on the same face
        if(const auto& shells = this->poly_con_.objects_on(fid))
        {
            // XXX: Here we assume that shells on a face is circular shells only
            for(const auto& sid : *shells)
            {
                if(filter(sid)) {continue;}

                auto sidp = this->get_shell(sid);
                const auto& sh = sidp.second.as_circular();
                // it is okay to use 3D distance because both are on the same face
                const Real dist = length(pos - sh.position()) - sh.shape().radius();
                if(dist < radius)
                {
                    retval.push_back(std::make_pair(std::move(sidp), dist));
                }
            }
        }

        // check shells on the neighborling faces
        for(const auto& nfid : polygon_->neighbor_faces_of(fid))
        {
            // XXX: Here we assume that shells on a face is circular shells only
            if(const auto& shells = this->poly_con_.objects_on(nfid))
            {
                for(const auto& sid : *shells)
                {
                    if(filter(sid)) {continue;}

                    auto sidp = this->get_shell(sid);
                    const auto& sh = sidp.second.as_circular();

                    const Real dist = ecell4::polygon::distance(*polygon_,
                        center, std::make_pair(sh.position(), sh.fid())
                        ) - sh.shape().radius();

                    if(dist < radius)
                    {
                        retval.push_back(std::make_pair(std::move(sidp), dist));
                    }
                }
            }
        }
        // check shells on the neighborling vertices
        for(const auto& vid : polygon_->neighbor_vertices_of(fid))
        {
            // XXX: Here we assume that shells on a vertex is conical shells only
            if(const auto& shells = this->poly_con_.objects_on(vid))
            {
                for(const auto& sid : *shells)
                {
                    if(filter(sid)) {continue;}

                    auto sidp = this->get_shell(sid);
                    const auto& sh = sidp.second.as_conical();

                    const Real dist = ecell4::polygon::distance(*polygon_,
                        center, std::make_pair(sh.position(), sh.vid())
                        ) - sh.shape().slant_height();

                    if(dist < radius)
                    {
                        retval.push_back(std::make_pair(std::move(sidp), dist));
                    }
                }
            }
        }
        std::sort(retval.begin(), retval.end(), utils::pair_second_element_comparator<
                  std::pair<ShellID, Shell>, Real>());
        return retval;
    }

private:

    rtree_type rtree_;
    polygon_container_type poly_con_;
    std::shared_ptr<Polygon> polygon_;
};

} // ngfrd
} // ecell4
#endif//ECELL4_NGFRD_SHELL_CONTAINER_HPP
