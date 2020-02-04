#ifndef ECELL4_SGFRD_SHELL_CONTAINER
#define ECELL4_SGFRD_SHELL_CONTAINER
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/Circle.hpp>
#include <ecell4/core/Cone.hpp>
#include <ecell4/sgfrd/StructureRegistrator.hpp>
#include <ecell4/sgfrd/distance_calculator.hpp>
#include <ecell4/sgfrd/Shell.hpp>

namespace ecell4
{
namespace sgfrd
{

class ShellContainer
{
public:

    typedef ::ecell4::Polygon polygon_type;
    typedef typename polygon_type::FaceID   FaceID;
    typedef typename polygon_type::EdgeID   EdgeID;
    typedef typename polygon_type::VertexID VertexID;

    typedef ecell4::Circle         circle_type;
    typedef ecell4::ConicalSurface conical_surface_type;
    typedef Shell<circle_type, FaceID>            circular_shell_type;
    typedef Shell<conical_surface_type, VertexID> conical_surface_shell_type;

    static const int circular_shell = 0;
    static const int conical_shell  = 1;
    typedef boost::variant<
                circular_shell_type,       // <- 0
                conical_surface_shell_type // <- 1
            > storage_type;
    typedef storage_type shell_type;
    typedef std::pair<ShellID, storage_type> shell_id_pair_type;
    typedef std::vector<shell_id_pair_type> container_type;
    typedef typename ecell4::utils::get_mapper_mf<ShellID, std::size_t>::type
        shell_id_to_index_map_type;

    typedef StructureRegistrator<ShellID, FaceID>   face_registrator_type;
    typedef StructureRegistrator<ShellID, VertexID> vertex_registrator_type;

    struct face_register_updater;
    struct vertex_register_updater;
    struct register_cleaner;

public:

    ShellContainer(const boost::shared_ptr<polygon_type>& poly)
        : polygon_(poly), face_registrator_(*poly), vertex_registrator_(*poly)
    {}
    ~ShellContainer(){}

    template<typename shellT>
    void add_shell(const ShellID& id, const shellT& sh, const FaceID&   fid);
    template<typename shellT>
    void add_shell(const ShellID& id, const shellT& sh, const VertexID& vid);

    // add a shell, with overlap checking.
    // Generally, it cannot be performed because Multi allows its shells to
    // overlap each other. These functions are for Single/Pair cases only.
    template<typename shellT>
    void check_add_shell(const ShellID& id, const shellT& sh, const FaceID&   fid,
                         const std::string& context);
    template<typename shellT>
    void check_add_shell(const ShellID& id, const shellT& sh, const VertexID& vid,
                         const std::string& context);

    template<typename shellT>
    void update_shell(const ShellID& id, const shellT& sh,
                      const FaceID& fid);
    template<typename shellT>
    void update_shell(const ShellID& id, const shellT& sh,
                      const VertexID& vid);

    storage_type const& get_shell(const ShellID& id) const noexcept
    {
        return container_.at(shell_id_to_index_map_.find(id)->second).second;
    }
    storage_type&       get_shell(const ShellID& id) noexcept
    {
        return container_.at(shell_id_to_index_map_[id]).second;
    }

    void remove_shell(const ShellID& shid);

    std::vector<ShellID> const& list_shells_on(const FaceID&   fid) const noexcept
    {
        return face_registrator_.elements_over(fid);
    }
    std::vector<ShellID> const& list_shells_on(const VertexID& vid) const noexcept
    {
        return vertex_registrator_.elements_over(vid);
    }

    std::size_t num_shells() const {return container_.size();}

    // calculate distance as 3D object
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real>>
    list_shells_within_radius(const Real3& pos, const Real radius) const
    {
        return list_shells_within_radius_impl(pos, radius,
                [](const ShellID&) noexcept -> bool {return false;});
    }
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real>>
    list_shells_within_radius(const Real3& pos, const Real radius,
                              const ShellID& ignore) const
    {
        return list_shells_within_radius_impl(pos, radius,
                [ignore](const ShellID& shid) noexcept -> bool {
                    return shid == ignore;
                });
    }
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real>>
    list_shells_within_radius(const Real3& pos, const Real radius,
                              const ShellID& ignore1,
                              const ShellID& ignore2) const
    {
        return this->list_shells_within_radius_impl(pos, radius,
                [ignore1, ignore2](const ShellID& shid) noexcept -> bool {
                    return shid == ignore1 || shid == ignore2;
                });
    }

    // calculate distance along the polygon.
    template<typename strID>
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real>>
    list_shells_within_radius(const std::pair<Real3, strID>& pos,
                              const Real radius) const
    {
        return this->list_shells_within_radius_impl(pos, radius,
                [](const ShellID&) noexcept -> bool {return false;});
    }
    template<typename strID>
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real>>
    list_shells_within_radius(const std::pair<Real3, strID>& pos,
                              const Real radius, const ShellID& ignore) const
    {
        return this->list_shells_within_radius_impl(pos, radius,
                [ignore](const ShellID& shid) noexcept -> bool {
                    return shid == ignore;
                });
    }
    template<typename strID>
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real>>
    list_shells_within_radius(const std::pair<Real3, strID>& pos,
                              const Real radius, const ShellID& ignore1,
                              const ShellID& ignore2) const
    {
        return this->list_shells_within_radius_impl(pos, radius,
                [ignore1, ignore2](const ShellID& shid) noexcept -> bool {
                    return shid == ignore1 || shid == ignore2;
                });
    }

    std::vector<shell_id_pair_type> list_shells() const {return container_;}

private:

    template<typename strID, typename F>
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real>>
    list_shells_within_radius_impl(const std::pair<Real3, strID>& pos,
                                   const Real radius, F ignore) const;

    // 3D version
    template<typename F>
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real>>
    list_shells_within_radius_impl(const Real3& pos,
                                   const Real radius, F ignore) const;

private:

    boost::shared_ptr<const polygon_type> polygon_;
    container_type             container_;
    shell_id_to_index_map_type shell_id_to_index_map_;
    face_registrator_type      face_registrator_;
    vertex_registrator_type    vertex_registrator_;
};

struct ShellContainer::register_cleaner
    : public boost::static_visitor<void>
{
    ShellContainer& scon;
    ShellID sid;

    register_cleaner(ShellContainer& self, const ShellID& si)
        : scon(self), sid(si){}

    template<typename shapeT>
    void operator()(const Shell<shapeT, FaceID>& sh) const
    {
        scon.face_registrator_.remove(sid);
        return;
    }

    template<typename shapeT>
    void operator()(const Shell<shapeT, VertexID>& sh) const
    {
        scon.vertex_registrator_.remove(sid);
        return;
    }
};

struct ShellContainer::face_register_updater
    : public boost::static_visitor<void>
{
    ShellContainer& scon;
    ShellID       sid;
    FaceID  fid;

    face_register_updater(ShellContainer& self, ShellID s, FaceID f)
        : scon(self), sid(s), fid(f)
    {}

    template<typename shapeT>
    void operator()(const Shell<shapeT, FaceID>& sh) const
    {
        scon.face_registrator_.update(sid, fid);
        return;
    }

    template<typename shapeT>
    void operator()(const Shell<shapeT, VertexID>& sh) const
    {
        scon.vertex_registrator_.remove(sid, sh.structure_id());
        scon.face_registrator_.emplace(sid, fid);
        return;
    }
};

struct ShellContainer::vertex_register_updater
    : public boost::static_visitor<void>
{
    ShellContainer& scon;
    ShellID         sid;
    VertexID  vid;

    vertex_register_updater(ShellContainer& self, ShellID s, VertexID v)
        : scon(self), sid(s), vid(v)
    {}

    template<typename shapeT>
    void operator()(const Shell<shapeT, FaceID>& sh) const
    {
        scon.face_registrator_.remove(sid, sh.structure_id());
        scon.vertex_registrator_.emplace(sid, vid);
        return;
    }

    template<typename shapeT>
    void operator()(const Shell<shapeT, VertexID>& sh) const
    {
        scon.vertex_registrator_.update(sid, vid);
        return;
    }
};

inline void ShellContainer::remove_shell(const ShellID& shid)
{
    if(shell_id_to_index_map_.count(shid) == 0)
    {
        throw std::invalid_argument("shellcontianer doesnt have the shell");
    }

    const std::size_t idx = shell_id_to_index_map_[shid];
    boost::apply_visitor(register_cleaner(*this, shid),
                         container_.at(idx).second);

    container_.at(idx) = container_.back();
    shell_id_to_index_map_[container_.back().first] = idx;
    container_.pop_back();
    return ;
}


template<typename shellT>
void ShellContainer::add_shell(
        const ShellID& id, const shellT& sh, const FaceID& fid)
{
    if(shell_id_to_index_map_.count(id) == 1)
    {
        throw std::invalid_argument("shellcontianer already have the shell");
    }

    const std::size_t idx = container_.size();
    shell_id_to_index_map_[id] = idx;
    face_registrator_.emplace(id, fid);
    container_.push_back(std::make_pair(id, storage_type(sh)));
    return;
}

template<typename shellT>
void ShellContainer::add_shell(
        const ShellID& id, const shellT& sh, const VertexID& vid)
{
    if(shell_id_to_index_map_.count(id) == 1)
    {
        throw std::invalid_argument("shellcontianer already have the shell");
    }
    const std::size_t idx = container_.size();
    shell_id_to_index_map_[id] = idx;
    vertex_registrator_.emplace(id, vid);
    container_.push_back(std::make_pair(id, storage_type(sh)));
    return;
}


template<typename shellT>
void ShellContainer::check_add_shell(
        const ShellID& id, const shellT& sh, const FaceID& fid,
        const std::string& context)
{
    if(shell_id_to_index_map_.count(id) == 1)
    {
        throw std::invalid_argument("shellcontianer already have the shell");
    }

    /* overlap check */ {
        std::vector<std::pair<std::pair<ShellID, storage_type>, Real>
            > ovlp = this->list_shells_within_radius(
                    std::make_pair(sh.position(), fid), sh.size());
        if(!ovlp.empty())
        {
            std::cerr << "WARNING: circular shells overlap!\n";
            std::cerr << "context: " << context << '\n';
            for(const auto& ov: ovlp)
            {
                std::cerr << "       : shell " << ov.first.first << " at "
                          << ov.second - sh.size() << "distant.\n";
            }
            std::cerr << std::flush;
        }
    }

    const std::size_t idx = container_.size();
    shell_id_to_index_map_[id] = idx;
    face_registrator_.emplace(id, fid);
    container_.push_back(std::make_pair(id, storage_type(sh)));
    return;
}

template<typename shellT>
void ShellContainer::check_add_shell(
        const ShellID& id, const shellT& sh, const VertexID& vid,
        const std::string& context)
{
    if(shell_id_to_index_map_.count(id) == 1)
    {
        throw std::invalid_argument("shellcontianer already have the shell");
    }

    /* overlap check */{
        std::vector<std::pair<std::pair<ShellID, storage_type>, Real>
            > ovlp = this->list_shells_within_radius(
                    std::make_pair(sh.position(), vid), sh.size());
        if(!ovlp.empty())
        {
            std::cerr << "WARNING: conical shells overlap!\n";
            std::cerr << "context: " << context << '\n';
            for(const auto& ov: ovlp)
            {
                std::cerr << "       : shell " << ov.first.first << " at "
                          << ov.second << "distant.\n";
            }
            std::cerr << std::flush;
        }
    }

    const std::size_t idx = container_.size();
    shell_id_to_index_map_[id] = idx;
    vertex_registrator_.emplace(id, vid);
    container_.push_back(std::make_pair(id, storage_type(sh)));
    return;
}


template<typename shellT>
void ShellContainer::update_shell(
        const ShellID& id, const shellT& sh, const FaceID& fid)
{
    if(shell_id_to_index_map_.count(id) == 0)
    {
        throw std::invalid_argument("shellcontianer doesnt have the shell");
    }

    const std::size_t idx = shell_id_to_index_map_[id];
    boost::apply_visitor(face_register_updater(*this, id, fid),
                         container_.at(idx).second);
    container_.at(idx).second = sh;
    return;
}

template<typename shellT>
void ShellContainer::update_shell(
        const ShellID& id, const shellT& sh, const VertexID& vid)
{
    if(shell_id_to_index_map_.count(id) == 0)
    {
        throw std::invalid_argument("shellcontianer doesnt have the shell");
    }

    const std::size_t idx = shell_id_to_index_map_[id];
    boost::apply_visitor(vertex_register_updater(*this, id, vid),
                         container_.at(idx).second);
    container_.at(idx).second = sh;
    return;
}

template<typename F>
std::vector<std::pair<std::pair<ShellID, ShellContainer::storage_type>, Real>>
ShellContainer::list_shells_within_radius_impl(
        const Real3& pos, const Real radius, F ignore) const
{
    throw NotSupported(
            "sGFRD does not support list_shells_within_radius in 3D space");
//     //XXX need more sophisticated way than brute-force searching
//     std::vector<std::pair<std::pair<ShellID, storage_type>, Real>> retval;
//     const distance_calculator distance(pos);
//     for(const auto& shp : this->container_)
//     {
//         const auto& shid  = shp.first;
//         const auto& shell = shp.second;
//         if(ignore(shid)) {continue;}
//
//         const Real dist = boost::apply_visitor(distance, shell);
//         if(dist < radius)
//         {
//             retval.emplace_back(shp, dist);
//         }
//     }
//     std::sort(retval.begin(), retval.end(),
//               ecell4::utils::pair_second_element_comparator<
//                   std::pair<ShellID, storage_type>, Real>());
//     return retval;
}

template<typename strID, typename F>
std::vector<std::pair<std::pair<ShellID, ShellContainer::storage_type>, Real>>
ShellContainer::list_shells_within_radius_impl(
        const std::pair<Real3, strID>& pos, const Real radius, F ignore) const
{
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real>> retval;
    const distance_calculator_on_surface<strID> distance_on_surf(pos, *polygon_);

    // check shells on the same position (either face or vertex)
    for(const ShellID& shid : list_shells_on(pos.second))
    {
        if(ignore(shid)) {continue;}

        const auto& shell = this->get_shell(shid);
        const Real  dist  = boost::apply_visitor(distance_on_surf, shell);
        if(dist < radius)
        {
            retval.emplace_back(std::make_pair(shid, shell), dist);
        }
    }

    const auto neighbor_faces = polygon_->neighbor_faces_of   (pos.second);
    const auto neighbor_vtxs  = polygon_->neighbor_vertices_of(pos.second);

    // check shells on the neighboring faces
    for(const FaceID& fid : neighbor_faces)
    {
        for(const ShellID& shid : list_shells_on(fid))
        {
            if(ignore(shid)) {continue;}

            const auto& shell = this->get_shell(shid);
            const Real  dist  = boost::apply_visitor(distance_on_surf, shell);
            if(dist < radius)
            {
                retval.emplace_back(std::make_pair(shid, shell), dist);
            }
        }
    }
    // check shells on the neighboring vertices
    for(const VertexID& vid : neighbor_vtxs)
    {
        for(const ShellID& shid : list_shells_on(vid))
        {
            if(ignore(shid)) {continue;}

            const auto& shell = this->get_shell(shid);
            const Real  dist  = boost::apply_visitor(distance_on_surf, shell);
            if(dist < radius)
            {
                retval.emplace_back(std::make_pair(shid, shell), dist);
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ShellID, storage_type>, Real>());

    // -----------------------------------------------------------------------
    // check double count.
    // neighbor_faces and neighbor_vtxs should not contain pos.second itself.
    // So if the polygon is okay, there will not be overlap in retval.

    std::set<ShellID> shellids;
    for(const auto& found : retval)
    {
        const auto& shid = found.first.first;
        if(shellids.count(shid) != 0)
        {
            std::ostringstream oss;
            oss << "Error: broken Polygon: Shell " << shid << " found twice.\n";
            oss << "neighboring faces of " << pos.second << " are";
            for(const FaceID& fid : neighbor_faces)
            {
                oss << ", " << fid;
            }
            oss << ".\n";
            oss << "neighboring vertices of " << pos.second << " are";
            for(const VertexID& vid : neighbor_vtxs)
            {
                oss << ", " << vid;
            }
            oss << ".\n";
            face_registrator_  .dump(std::cerr);
            vertex_registrator_.dump(std::cerr);
            throw std::runtime_error(oss.str());
        }
        shellids.insert(found.first.first);
    }
    return retval;
}

}// sgfrd
}// ecell4
#endif// ECELL4_SGFRD_SHELL_CONTAINER
