#ifndef ECELL4_SGFRD_SHELL_CONTAINER
#define ECELL4_SGFRD_SHELL_CONTAINER
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/Circle.hpp>
#include <ecell4/core/Cone.hpp>
#include <ecell4/sgfrd/Shell.hpp>

namespace ecell4
{
namespace sgfrd
{

template<typename T_polygon_traits>
class ShellContainer
{
public:
    typedef T_polygon_traits traits_type;
    typedef Polygon<traits_type> polygon_type;
    typedef typename polygon_type::face_id_type   face_id_type;
    typedef typename polygon_type::edge_id_type   edge_id_type;
    typedef typename polygon_type::vertex_id_type vertex_id_type;

    typedef ecell4::Circle         circle_type;
    typedef ecell4::ConicalSurface conical_surface_type;
    typedef Shell<circle_type, face_id_type>            circular_shell_type;
    typedef Shell<conical_surface_type, vertex_id_type> conical_surface_shell_type;

    typedef boost::variant<
                circular_shell_type,
                conical_surface_shell_type
            > storage_type;
    typedef std::pair<ShellID, storage_type> shell_id_pair_type;
    typedef std::vector<shell_id_pair_type> container_type;
    typedef typename ecell4::utils::get_mapper_mf<ShellID, std::size_t>::type
        shell_id_to_index_map_type;

    typedef StructureRegistrator<ShellID, face_id_type, traits_type>
        face_registrator_type;
    typedef StructureRegistrator<ShellID, vertex_id_type, traits_type>
        vertex_registrator_type;

    template<typename sidT>
    struct register_updater;
    struct register_cleaner;

public:

    ShellContainer(const polygon_type& poly)
        : polygon_(poly), face_registrator_(poly), vertex_registrator_(poly)
    {}
    ~ShellContainer(){}

    template<typename shellT>
    void add_shell(const ShellID& id, const shellT& sh, const face_id_type&   fid);
    template<typename shellT>
    void add_shell(const ShellID& id, const shellT& sh, const vertex_id_type& vid);

    template<typename shellT>
    void update_shell(const ShellID& id, const shellT& sh,
                      const face_id_type& fid) const;
    void update_shell(const ShellID& id, const shellT& sh,
                      const vertex_id_type& vid) const;

    storage_type const& get_shell(const ShellID& id) const;
    storage_type&       get_shell(const ShellID& id);

    void remove_shell(const ShellID& id);

    // calculate distance as 3D object
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const Real3& pos, const Real radius);
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const Real3& pos, const Real radius,
            const ShellID& ignore);
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const Real3& pos, const Real radius,
            const ShellID& ignore1, const ShellID& ignore2);

    //calculate distance between position on face
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const std::pair<Real3, face_id_type>& pos,
            const Real radius);
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const std::pair<Real3, face_id_type>& pos,
            const Real radius, const ShellID& ignore);
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const std::pair<Real3, face_id_type>& pos,
            const Real radius, const ShellID& ignore1, const ShellID& ignore2);

    //calculate distance between position on vertex
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const std::pair<Real3, vertex_id_type>& pos,
            const Real radius);
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const std::pair<Real3, vertex_id_type>& pos,
            const Real radius, const ShellID& ignore);
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const std::pair<Real3, vertex_id_type>& pos,
            const Real radius, const ShellID& ignore1, const ShellID& ignore2);

private:

    const polygon_type&        polygon_;
    container_type             container_;
    shell_id_to_index_map_type shell_id_to_index_map_;
    face_registrator_type      face_registrator_;
    vertex_registrator_type    vertex_registrator_;
};

template<typename T_pt>
struct ShellContainer<T_pt>::register_cleaner : public boost::static_visitor<void>
{
    ShellContainer<T_pt>& scon;
    register_cleaner(ShellContainer<T_pt>& self) : scon(self){}

    template<typename shapeT>
    void operator()(const Shell<shapeT, face_id_type>& sh) const
    {
        scon.face_registrator_.remove(sh.structure_id());
        return;
    }

    template<typename shapeT>
    void operator()(const Shell<shapeT, vertex_id_type>& sh) const
    {
        scon.vertex_registrator_.remove(sh.structure_id());
        return;
    }
};

template<typename T_pt>
struct ShellContainer<T_pt>::register_updater<face_id_type>
    : public boost::static_visitor<void>
{
    ShellContainer<T_pt>& scon;
    ShellID      sid;
    face_id_type fid;

    register_cleaner(ShellContainer<T_pt>& self, ShellID s, face_id_type f)
        : scon(self), sid(s), fid(f)
    {}

    template<typename shapeT>
    void operator()(const Shell<shapeT, face_id_type>& sh) const
    {
        scon.face_registrator_.update(sid, fid);
        return;
    }

    template<typename shapeT>
    void operator()(const Shell<shapeT, vertex_id_type>& sh) const
    {
        scon.vertex_registrator_.remove(sh.structure_id());
        scon.face_registrator_.emplace(sid, fid);
        return;
    }
};

template<typename T_pt>
struct ShellContainer<T_pt>::register_updater<vertex_id_type>
    : public boost::static_visitor<void>
{
    ShellContainer<T_pt>& scon;
    ShellID        sid;
    vertex_id_type vid;

    register_cleaner(ShellContainer<T_pt>& self, ShellID s, vertex_id_type v)
        : scon(self), sid(s), vid(v)
    {}

    template<typename shapeT>
    void operator()(const Shell<shapeT, face_id_type>& sh) const
    {
        scon.face_registrator_.remove(sh.structure_id());
        scon.vertex_registrator_.emplace(sid, vid);
        return;
    }

    template<typename shapeT>
    void operator()(const Shell<shapeT, vertex_id_type>& sh) const
    {
        scon.vertex_registrator_.update(sid, vid);
        return;
    }
};

template<typename T_pt>
template<typename shellT>
void ShellContainer<T_pt>::add_shell(
        const ShellID& id, const shellT& sh, const face_id_type& fid)
{
    if(shell_id_to_index_map_.count(id) == 1)
        throw std::invalid_argument("shellcontianer already have the shell");
    const std::size_t idx = container_.size();
    shell_id_to_index_map_[id] = idx;
    face_registrator_.emplace(id, fid);
    container_.push_back(std::make_pair(id, storage_type(sh)));
    return;
}

template<typename T_pt>
template<typename shellT>
void ShellContainer<T_pt>::add_shell(
        const ShellID& id, const shellT& sh, const vertex_id_type& vid)
{
    if(shell_id_to_index_map_.count(id) == 1)
        throw std::invalid_argument("shellcontianer already have the shell");
    const std::size_t idx = container_.size();
    shell_id_to_index_map_[id] = idx;
    vertex_registrator_.emplace(id, vid);
    container_.push_back(std::make_pair(id, storage_type(sh)));
    return;
}

template<typename T_pt>
template<typename shellT>
void ShellContainer<T_pt>::update_shell(
        const ShellID& id, const shellT& sh, const face_id_type& fid)
{
    if(shell_id_to_index_map.count(id) == 0)
        throw std::invalid_argument("shellcontianer doesnt have the shell");
    const std::size_t idx = shell_id_to_index_map_[id];
    boost::apply_visitor(register_updator<face_id_type>(*this, id, fid),
                         container_.at(idx));
    container_.at(idx).second = sh;
    return;
}

template<typename T_pt>
template<typename shellT>
void ShellContainer<T_pt>::update_shell(
        const ShellID& id, const shellT& sh, const vertex_id_type& vid)
{
    if(shell_id_to_index_map.count(id) == 0)
        throw std::invalid_argument("shellcontianer doesnt have the shell");
    const std::size_t idx = shell_id_to_index_map_[id];
    boost::apply_visitor(register_updator<vertex_id_type>(*this, id, vid),
                         container_.at(idx));
    container_.at(idx).second = sh;
    return;
}

template<typename T_pt>
typename ShellContainer<T_pt>::storage_type const&
ShellContainer<T_pt>::get_shell(const ShellID& id) const
{
    return container_.at(shell_id_to_index_map_.find(id)->second).second;
}

template<typename T_pt>
typename ShellContainer<T_pt>::storage_type&
ShellContainer<T_pt>::get_shell(const ShellID& id)
{
    return container_.at(shell_id_to_index_map_[id]).second;
}

template<typename T_pt>
void ShellContainer<T_pt>::remove_shell(const ShellID& id)
{
    if(shell_id_to_index_map.count(id) == 0)
        throw std::invalid_argument("shellcontianer doesnt have the shell");
    const std::size_t idx = shell_id_to_index_map_[id];
    boost::apply_visitor(register_cleaner(*this), container_.at(idx));

    container_.at(idx) = container_.back();
    container_.pop_back();
    return ;
}


}// sgfrd
}// ecell4
#endif// ECELL4_SGFRD_SHELL_CONTAINER
