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
    typedef storage_type shell_type;
    typedef std::pair<ShellID, storage_type> shell_id_pair_type;
    typedef std::vector<shell_id_pair_type> container_type;
    typedef typename ecell4::utils::get_mapper_mf<ShellID, std::size_t>::type
        shell_id_to_index_map_type;

    typedef StructureRegistrator<ShellID, face_id_type, traits_type>
        face_registrator_type;
    typedef StructureRegistrator<ShellID, vertex_id_type, traits_type>
        vertex_registrator_type;

    struct face_register_updater;
    struct vertex_register_updater;
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
                      const face_id_type& fid);
    template<typename shellT>
    void update_shell(const ShellID& id, const shellT& sh,
                      const vertex_id_type& vid);

    storage_type const& get_shell(const ShellID& id) const;
    storage_type&       get_shell(const ShellID& id);

    void remove_shell(const ShellID& id);

    std::vector<ShellID> const& list_shells_on(const face_id_type&) const;
    std::vector<ShellID> const& list_shells_on(const vertex_id_type&) const;

    std::size_t num_shells() const {return container_.size();}

    // calculate distance as 3D object
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const Real3& pos, const Real radius) const;
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const Real3& pos, const Real radius,
            const ShellID& ignore) const;
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const Real3& pos, const Real radius,
            const ShellID& ignore1, const ShellID& ignore2) const;

    //calculate distance along the polygon
    template<typename strID>
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const std::pair<Real3, strID>& pos,
            const Real radius) const;
    template<typename strID>
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const std::pair<Real3, strID>& pos,
            const Real radius, const ShellID& ignore) const;
    template<typename strID>
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> >
        list_shells_within_radius(const std::pair<Real3, strID>& pos,
            const Real radius, const ShellID& ignore1, const ShellID& ignore2) const;

private:

    const polygon_type&        polygon_;
    container_type             container_;
    shell_id_to_index_map_type shell_id_to_index_map_;
    face_registrator_type      face_registrator_;
    vertex_registrator_type    vertex_registrator_;
};

template<typename T_pt>
struct ShellContainer<T_pt>::register_cleaner
    : public boost::static_visitor<void>
{
    ShellContainer<T_pt>& scon;
    ShellID sid;

    register_cleaner(ShellContainer<T_pt>& self, const ShellID& si)
        : scon(self), sid(si){}

    template<typename shapeT>
    void operator()(const Shell<shapeT, face_id_type>& sh) const
    {
        scon.face_registrator_.remove(sid);
        return;
    }

    template<typename shapeT>
    void operator()(const Shell<shapeT, vertex_id_type>& sh) const
    {
        scon.vertex_registrator_.remove(sid);
        return;
    }
};

template<typename T_pt>
struct ShellContainer<T_pt>::face_register_updater
    : public boost::static_visitor<void>
{
    ShellContainer<T_pt>& scon;
    ShellID      sid;
    face_id_type fid;

    face_register_updater(ShellContainer<T_pt>& self, ShellID s, face_id_type f)
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
struct ShellContainer<T_pt>::vertex_register_updater
    : public boost::static_visitor<void>
{
    ShellContainer<T_pt>& scon;
    ShellID        sid;
    vertex_id_type vid;

    vertex_register_updater(ShellContainer<T_pt>& self, ShellID s, vertex_id_type v)
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
    if(shell_id_to_index_map_.count(id) == 0)
        throw std::invalid_argument("shellcontianer doesnt have the shell");
    const std::size_t idx = shell_id_to_index_map_[id];
    boost::apply_visitor(face_register_updater(*this, id, fid),
                         container_.at(idx));
    container_.at(idx).second = sh;
    return;
}

template<typename T_pt>
template<typename shellT>
void ShellContainer<T_pt>::update_shell(
        const ShellID& id, const shellT& sh, const vertex_id_type& vid)
{
    if(shell_id_to_index_map_.count(id) == 0)
        throw std::invalid_argument("shellcontianer doesnt have the shell");
    const std::size_t idx = shell_id_to_index_map_[id];
    boost::apply_visitor(vertex_register_updater(*this, id, vid),
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
    if(shell_id_to_index_map_.count(id) == 0)
        throw std::invalid_argument("shellcontianer doesnt have the shell");
    const std::size_t idx = shell_id_to_index_map_[id];
    boost::apply_visitor(register_cleaner(*this, id),
                         container_.at(idx).second);

    container_.at(idx) = container_.back();
    shell_id_to_index_map_[container_.back().first] = idx;
    container_.pop_back();
    return ;
}

template<typename T_pt>
inline std::vector<ShellID> const&
ShellContainer<T_pt>::list_shells_on(const face_id_type& fid) const
{
    return face_registrator_.elements_over(fid);
}

template<typename T_pt>
inline std::vector<ShellID> const&
ShellContainer<T_pt>::list_shells_on(const vertex_id_type& vid) const
{
    return vertex_registrator_.elements_over(vid);
}

template<typename T_pt>
std::vector<std::pair<
    std::pair<ShellID, typename ShellContainer<T_pt>::storage_type>, Real> >
ShellContainer<T_pt>::list_shells_within_radius(
        const Real3& pos, const Real radius) const
{
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> > retval;
    const distance_calculator distance(pos);
    //XXX need more sophisticated way than brute-force searching

    for(typename container_type::const_iterator
        iter = container_.begin(); iter != container_.end(); ++iter)
    {
        const Real dist = boost::apply_visitor(distance, iter->second);
        if(dist < radius)
        {
            retval.push_back(std::make_pair(*iter, dist));
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ShellID, storage_type>, Real>());
    return retval;
}

template<typename T_pt>
std::vector<std::pair<
    std::pair<ShellID, typename ShellContainer<T_pt>::storage_type>, Real> >
ShellContainer<T_pt>::list_shells_within_radius(
        const Real3& pos, const Real radius, const ShellID& ignore) const
{
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> > retval;
    const distance_calculator distance(pos);
    //XXX need more sophisticated way than brute-force searching

    for(typename container_type::const_iterator
        iter = container_.begin(); iter != container_.end(); ++iter)
    {
        if(iter->first == ignore) continue;
        const Real dist = boost::apply_visitor(distance, iter->second);
        if(dist < radius)
        {
            retval.push_back(std::make_pair(*iter, dist));
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ShellID, storage_type>, Real>());
    return retval;
}

template<typename T_pt>
std::vector<std::pair<
    std::pair<ShellID, typename ShellContainer<T_pt>::storage_type>, Real> >
ShellContainer<T_pt>::list_shells_within_radius(
        const Real3& pos, const Real radius,
        const ShellID& ignore1, const ShellID& ignore2) const
{
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> > retval;
    const distance_calculator distance(pos);
    //XXX need more sophisticated way than brute-force searching

    for(typename container_type::const_iterator
        iter = container_.begin(); iter != container_.end(); ++iter)
    {
        if(iter->first == ignore1 || iter->first == ignore2) continue;
        const Real dist = boost::apply_visitor(distance, iter->second);
        if(dist < radius)
        {
            retval.push_back(std::make_pair(*iter, dist));
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ShellID, storage_type>, Real>());
    return retval;
}

template<typename T_pt>
template<typename strID>
std::vector<std::pair<
    std::pair<ShellID, typename ShellContainer<T_pt>::storage_type>, Real> >
ShellContainer<T_pt>::list_shells_within_radius(
        const std::pair<Real3, strID>& pos, const Real radius) const
{
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> > retval;
    const distance_calculator_on_surface<T_pt, strID>
        distance(pos, this->polygon_);

    std::vector<face_id_type>   neighborf = polygon_.at(pos.second).neighbor_faces;
    std::vector<vertex_id_type> neighborv = polygon_.at(pos.second).neighbor_vertices;

    for(typename std::vector<face_id_type>::const_iterator
        iter = neighborf.begin(); iter != neighborf.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            iter = shells.begin(); iter != shells.end(); ++iter)
        {
            const storage_type shell(this->get_shell(*iter));
            const Real dist = boost::apply_visitor(distance, shell);
            if(dist < radius)
            {
                retval.push_back(std::make_pair(
                            std::make_pair(*iter, shell), dist));
            }
        }
    }

    for(typename std::vector<vertex_id_type>::const_iterator
        iter = neighborv.begin(); iter != neighborv.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            iter = shells.begin(); iter != shells.end(); ++iter)
        {
            const storage_type shell(this->get_shell(*iter));
            const Real dist = boost::apply_visitor(distance, this->get_shell(*iter));
            if(dist < radius)
            {
                retval.push_back(std::make_pair(
                            std::make_pair(*iter, shell), dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ShellID, storage_type>, Real>());
    return retval;
}

template<typename T_pt>
template<typename strID>
std::vector<std::pair<
    std::pair<ShellID, typename ShellContainer<T_pt>::storage_type>, Real> >
ShellContainer<T_pt>::list_shells_within_radius(
        const std::pair<Real3, strID>& pos, const Real radius,
        const ShellID& ignore) const
{
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> > retval;
    const distance_calculator_on_surface<T_pt, strID>
        distance(pos, this->polygon_);

    std::vector<face_id_type>   neighborf = polygon_.at(pos.second).neighbor_faces;
    std::vector<vertex_id_type> neighborv = polygon_.at(pos.second).neighbor_vertices;

    for(typename std::vector<face_id_type>::const_iterator
        iter = neighborf.begin(); iter != neighborf.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            iter = shells.begin(); iter != shells.end(); ++iter)
        {
            if(*iter == ignore) continue;
            const storage_type shell(this->get_shell(*iter));
            const Real dist = boost::apply_visitor(distance, shell);
            if(dist < radius)
            {
                retval.push_back(std::make_pair(
                            std::make_pair(*iter, shell), dist));
            }
        }
    }

    for(typename std::vector<vertex_id_type>::const_iterator
        iter = neighborv.begin(); iter != neighborv.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            iter = shells.begin(); iter != shells.end(); ++iter)
        {
            if(*iter == ignore) continue;
            const storage_type shell(this->get_shell(*iter));
            const Real dist = boost::apply_visitor(distance, shell);
            if(dist < radius)
            {
                retval.push_back(std::make_pair(
                            std::make_pair(*iter,shell), dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ShellID, storage_type>, Real>());
    return retval;
}

template<typename T_pt>
template<typename strID>
std::vector<std::pair<
    std::pair<ShellID, typename ShellContainer<T_pt>::storage_type>, Real> >
ShellContainer<T_pt>::list_shells_within_radius(
        const std::pair<Real3, strID>& pos, const Real radius,
        const ShellID& ignore1, const ShellID& ignore2) const
{
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> > retval;
    const distance_calculator_on_surface<T_pt, strID>
        distance(pos, this->polygon_);

    std::vector<face_id_type>   neighborf = polygon_.at(pos.second).neighbor_faces;
    std::vector<vertex_id_type> neighborv = polygon_.at(pos.second).neighbor_vertices;

    for(typename std::vector<face_id_type>::const_iterator
        iter = neighborf.begin(); iter != neighborf.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            iter = shells.begin(); iter != shells.end(); ++iter)
        {
            if(*iter == ignore1 || *iter == ignore2) continue;
            const storage_type shell(this->get_shell(*iter));
            const Real dist = boost::apply_visitor(distance, shell);
            if(dist < radius)
            {
                retval.push_back(std::make_pair(
                            std::make_pair(*iter, shell), dist));
            }
        }
    }

    for(typename std::vector<vertex_id_type>::const_iterator
        iter = neighborv.begin(); iter != neighborv.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            iter = shells.begin(); iter != shells.end(); ++iter)
        {
            if(*iter == ignore1 || *iter == ignore2) continue;

            const storage_type shell(this->get_shell(*iter));
            const Real dist = boost::apply_visitor(distance, shell);
            if(dist < radius)
            {
                retval.push_back(std::make_pair(
                            std::make_pair(*iter, shell), dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ShellID, storage_type>, Real>());
    return retval;
}


}// sgfrd
}// ecell4
#endif// ECELL4_SGFRD_SHELL_CONTAINER
