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

    typedef StructureRegistrator<ShellID, FaceID, traits_type>
        face_registrator_type;
    typedef StructureRegistrator<ShellID, VertexID, traits_type>
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
    void add_shell(const ShellID& id, const shellT& sh, const FaceID&   fid);
    template<typename shellT>
    void add_shell(const ShellID& id, const shellT& sh, const VertexID& vid);

    template<typename shellT>
    void update_shell(const ShellID& id, const shellT& sh,
                      const FaceID& fid);
    template<typename shellT>
    void update_shell(const ShellID& id, const shellT& sh,
                      const VertexID& vid);

    storage_type const& get_shell(const ShellID& id) const;
    storage_type&       get_shell(const ShellID& id);

    void remove_shell(const ShellID& id);

    std::vector<ShellID> const& list_shells_on(const FaceID&) const;
    std::vector<ShellID> const& list_shells_on(const VertexID&) const;

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

    std::vector<shell_id_pair_type>
    list_shells() const {return container_;}

private:

    const polygon_type&        polygon_;
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

template<typename shellT>
void ShellContainer::add_shell(
        const ShellID& id, const shellT& sh, const FaceID& fid)
{
    if(shell_id_to_index_map_.count(id) == 1)
        throw std::invalid_argument("shellcontianer already have the shell");

//     /* overlap check */{
//         std::vector<std::pair<std::pair<ShellID, storage_type>, Real>
//             > ovlp = this->list_shells_within_radius(
//                     std::make_pair(sh.position(), fid), sh.size());
//         if(!ovlp.empty())
//         {
//             std::cout << "WARNING: circular shells overlap!" << std::endl;
//         }
//     }

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
        throw std::invalid_argument("shellcontianer already have the shell");

//     /* overlap check */{
//         std::vector<std::pair<std::pair<ShellID, storage_type>, Real>
//             > ovlp = this->list_shells_within_radius(
//                     std::make_pair(sh.position(), vid), sh.size());
//         if(!ovlp.empty())
//         {
//             std::cout << "WARNING: conical shells overlap!" << std::endl;
//         }
//     }

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
        throw std::invalid_argument("shellcontianer doesnt have the shell");

//     /* overlap check */{
//         std::vector<std::pair<std::pair<ShellID, storage_type>, Real>
//             > ovlp = this->list_shells_within_radius(
//                     std::make_pair(sh.position(), fid), sh.size(), id);
//         if(!ovlp.empty())
//         {
//             std::cout << "WARNING: circular shells overlap!" << std::endl;
//         }
//     }

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
        throw std::invalid_argument("shellcontianer doesnt have the shell");

//     /* overlap check */{
//         std::vector<std::pair<std::pair<ShellID, storage_type>, Real>
//             > ovlp = this->list_shells_within_radius(
//                     std::make_pair(sh.position(), vid), sh.size(), id);
//         if(!ovlp.empty())
//         {
//             std::cout << "WARNING: circular shell overlaps!" << std::endl;
//         }
//     }

    const std::size_t idx = shell_id_to_index_map_[id];
    boost::apply_visitor(vertex_register_updater(*this, id, vid),
                         container_.at(idx).second);
    container_.at(idx).second = sh;
    return;
}

inline ShellContainer::storage_type const&
ShellContainer::get_shell(const ShellID& id) const
{
    return container_.at(shell_id_to_index_map_.find(id)->second).second;
}

inline ShellContainer::storage_type&
ShellContainer::get_shell(const ShellID& id)
{
    return container_.at(shell_id_to_index_map_[id]).second;
}

inline void ShellContainer::remove_shell(const ShellID& id)
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

inline std::vector<ShellID> const&
ShellContainer::list_shells_on(const FaceID& fid) const
{
    return face_registrator_.elements_over(fid);
}

inline std::vector<ShellID> const&
ShellContainer::list_shells_on(const VertexID& vid) const
{
    return vertex_registrator_.elements_over(vid);
}

inline std::vector<std::pair<
    std::pair<ShellID, ShellContainer::storage_type>, Real> >
ShellContainer::list_shells_within_radius(
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

inline std::vector<std::pair<
    std::pair<ShellID, ShellContainer::storage_type>, Real> >
ShellContainer::list_shells_within_radius(
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

std::vector<std::pair<
    std::pair<ShellID, ShellContainer::storage_type>, Real> >
ShellContainer::list_shells_within_radius(
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

template<typename strID>
std::vector<std::pair<
    std::pair<ShellID, ShellContainer::storage_type>, Real> >
ShellContainer::list_shells_within_radius(
        const std::pair<Real3, strID>& pos, const Real radius) const
{
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> > retval;
    const distance_calculator_on_surface<T_pt, strID>
        distance_on_surf(pos, this->polygon_);

    {
        const std::vector<ShellID>& shells = list_shells_on(pos.second);
        for(typename std::vector<ShellID>::const_iterator
            jter = shells.begin(); jter != shells.end(); ++jter)
        {
            const storage_type shell(this->get_shell(*jter));
            const Real dist = boost::apply_visitor(distance_on_surf, shell);
            if(dist < radius)
            {
                retval.push_back(std::make_pair(
                            std::make_pair(*jter, shell), dist));
            }
        }
    }

    const std::vector<FaceID>   neighborf = polygon_.at(pos.second).neighbor_faces;
    const std::vector<VertexID> neighborv = polygon_.at(pos.second).neighbor_vertices;

    for(typename std::vector<FaceID>::const_iterator
        iter = neighborf.begin(); iter != neighborf.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            jter = shells.begin(); jter != shells.end(); ++jter)
        {
            const storage_type shell(this->get_shell(*jter));
            const Real dist = boost::apply_visitor(distance_on_surf, shell);
            if(dist < radius)
            {
                for(std::size_t i=0; i<retval.size(); ++i)
                {
                    // to avoid double-count;
                    // to check the consistency of polygon and structure_registrator.
                    if(retval.at(i).first.first == *jter)
                    {
                        std::cerr << "Error: broken consistency in structure id\n";
                        std::cerr << "ShellContainer::list_shells_within_radius\n";
                        std::cerr << "FaceID " << pos.second << " has ...\n";
                        for(typename std::vector<FaceID>::const_iterator
                            nfi(neighborf.begin()), nfe(neighborf.end());
                            nfi != nfe; ++nfi)
                        {
                            std::cerr << "  FaceID = " << *nfi << std::endl;
                        }
                        std::cerr << "neighbors.\n";
                        std::cerr << "Shell " << *jter
                                  << " is found on both Face " << *iter
                                  << " and Face " << pos.second << std::endl;
                        face_registrator_.dump(std::cerr);
                        assert(false);
                    }
                }

                retval.push_back(std::make_pair(
                            std::make_pair(*jter, shell), dist));
            }
        }
    }

    for(typename std::vector<VertexID>::const_iterator
        iter = neighborv.begin(); iter != neighborv.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            jter = shells.begin(); jter != shells.end(); ++jter)
        {
            const storage_type shell(this->get_shell(*jter));
            const Real dist = boost::apply_visitor(distance_on_surf, this->get_shell(*jter));
            if(dist < radius)
            {
                for(std::size_t i=0; i<retval.size(); ++i)
                {
                    if(retval.at(i).first.first == *jter)
                    {
                        std::cerr << "Error: broken consistency in structure id\n";
                        std::cerr << "ShellContainer::list_shells_within_radius\n";
                        std::cerr << "Vertex " << pos.second << " has ...\n";
                        for(typename std::vector<VertexID>::const_iterator
                            nvi(neighborv.begin()), nve(neighborv.end());
                            nvi != nve; ++nvi)
                        {
                            std::cerr << "  VertexID = " << *nvi << std::endl;
                        }
                        std::cerr << "neighbors.\n";
                        std::cerr << "Shell " << *jter
                                  << " is found on both Vertex " << *iter
                                  << " and another Face." << std::endl;
                        face_registrator_.dump(std::cerr);
                        vertex_registrator_.dump(std::cerr);
                        assert(false);
                    }
                }
                retval.push_back(std::make_pair(
                            std::make_pair(*jter, shell), dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ShellID, storage_type>, Real>());
    return retval;
}

template<typename strID>
std::vector<std::pair<
    std::pair<ShellID, ShellContainer::storage_type>, Real> >
ShellContainer::list_shells_within_radius(
        const std::pair<Real3, strID>& pos, const Real radius,
        const ShellID& ignore) const
{
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> > retval;
    const distance_calculator_on_surface<T_pt, strID>
        distance_on_surf(pos, this->polygon_);

    {
        const std::vector<ShellID>& shells = list_shells_on(pos.second);
        for(typename std::vector<ShellID>::const_iterator
            jter = shells.begin(); jter != shells.end(); ++jter)
        {
            const storage_type shell(this->get_shell(*jter));
            const Real dist = boost::apply_visitor(distance_on_surf, shell);
            if(dist < radius)
            {
                retval.push_back(std::make_pair(
                            std::make_pair(*jter, shell), dist));
            }
        }
    }

    std::vector<FaceID>   neighborf = polygon_.at(pos.second).neighbor_faces;
    std::vector<VertexID> neighborv = polygon_.at(pos.second).neighbor_vertices;

    for(typename std::vector<FaceID>::const_iterator
        iter = neighborf.begin(); iter != neighborf.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            jter = shells.begin(); jter != shells.end(); ++jter)
        {
            if(*jter == ignore) continue;
            const storage_type shell(this->get_shell(*jter));
            const Real dist = boost::apply_visitor(distance_on_surf, shell);
            if(dist < radius)
            {
                for(std::size_t i=0; i<retval.size(); ++i)
                {
                    if(retval.at(i).first.first == *jter)
                    {
                        std::cerr << "Error: broken consistency in structure id\n";
                        std::cerr << "ShellContainer::list_shells_within_radius\n";
                        std::cerr << "FaceID " << pos.second << " has ...\n";
                        for(typename std::vector<FaceID>::const_iterator
                            nfi(neighborf.begin()), nfe(neighborf.end());
                            nfi != nfe; ++nfi)
                        {
                            std::cerr << "  FaceID = " << *nfi << std::endl;
                        }
                        std::cerr << "neighbors.\n";
                        std::cerr << "Shell " << *jter
                                  << " is found on both Face " << *iter
                                  << " and Face " << pos.second << std::endl;
                        face_registrator_.dump(std::cerr);
                        assert(false);
                    }
                }

                retval.push_back(std::make_pair(
                            std::make_pair(*jter, shell), dist));
            }
        }
    }

    for(typename std::vector<VertexID>::const_iterator
        iter = neighborv.begin(); iter != neighborv.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            jter = shells.begin(); jter != shells.end(); ++jter)
        {
            if(*jter == ignore) continue;
            const storage_type shell(this->get_shell(*jter));
            const Real dist = boost::apply_visitor(distance_on_surf, shell);
            if(dist < radius)
            {
                for(std::size_t i=0; i<retval.size(); ++i)
                {
                    if(retval.at(i).first.first == *jter)
                    {
                        std::cerr << "Error: broken consistency in structure id\n";
                        std::cerr << "ShellContainer::list_shells_within_radius\n";
                        std::cerr << "Vertex " << pos.second << " has ...\n";
                        for(typename std::vector<VertexID>::const_iterator
                            nvi(neighborv.begin()), nve(neighborv.end());
                            nvi != nve; ++nvi)
                        {
                            std::cerr << "  VertexID = " << *nvi << std::endl;
                        }
                        std::cerr << "neighbors.\n";
                        std::cerr << "Shell " << *jter
                                  << " is found on both Vertex " << *iter
                                  << " and another Face." << std::endl;
                        face_registrator_.dump(std::cerr);
                        vertex_registrator_.dump(std::cerr);
                        assert(false);
                    }
                }

                retval.push_back(std::make_pair(
                            std::make_pair(*jter,shell), dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ShellID, storage_type>, Real>());
    return retval;
}

template<typename strID>
std::vector<std::pair<
    std::pair<ShellID, ShellContainer::storage_type>, Real> >
ShellContainer::list_shells_within_radius(
        const std::pair<Real3, strID>& pos, const Real radius,
        const ShellID& ignore1, const ShellID& ignore2) const
{
    std::vector<std::pair<std::pair<ShellID, storage_type>, Real> > retval;
    const distance_calculator_on_surface<T_pt, strID>
        distance_on_surf(pos, this->polygon_);

    {
        const std::vector<ShellID>& shells = list_shells_on(pos.second);
        for(typename std::vector<ShellID>::const_iterator
            jter = shells.begin(); jter != shells.end(); ++jter)
        {
            const storage_type shell(this->get_shell(*jter));
            const Real dist = boost::apply_visitor(distance_on_surf, shell);
            if(dist < radius)
            {
                retval.push_back(std::make_pair(
                            std::make_pair(*jter, shell), dist));
            }
        }
    }

    std::vector<FaceID>   neighborf = polygon_.at(pos.second).neighbor_faces;
    std::vector<VertexID> neighborv = polygon_.at(pos.second).neighbor_vertices;

    for(typename std::vector<FaceID>::const_iterator
        iter = neighborf.begin(); iter != neighborf.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            jter = shells.begin(); jter != shells.end(); ++jter)
        {
            if(*jter == ignore1 || *jter == ignore2) continue;
            const storage_type shell(this->get_shell(*jter));
            const Real dist = boost::apply_visitor(distance_on_surf, shell);
            if(dist < radius)
            {
                for(std::size_t i=0; i<retval.size(); ++i)
                {
                    if(retval.at(i).first.first == *jter)
                    {
                        std::cerr << "Error: broken consistency in structure id\n";
                        std::cerr << "ShellContainer::list_shells_within_radius\n";
                        std::cerr << "FaceID " << pos.second << " has ...\n";
                        for(typename std::vector<FaceID>::const_iterator
                            nfi(neighborf.begin()), nfe(neighborf.end());
                            nfi != nfe; ++nfi)
                        {
                            std::cerr << "  FaceID = " << *nfi << std::endl;
                        }
                        std::cerr << "neighbors.\n";
                        std::cerr << "Shell " << *jter
                                  << " is found on both Face " << *iter
                                  << " and Face " << pos.second << std::endl;
                        face_registrator_.dump(std::cerr);
                        assert(false);
                    }
                }

                retval.push_back(std::make_pair(
                            std::make_pair(*jter, shell), dist));
            }
        }
    }

    for(typename std::vector<VertexID>::const_iterator
        iter = neighborv.begin(); iter != neighborv.end(); ++iter)
    {
        const std::vector<ShellID>& shells = list_shells_on(*iter);
        for(typename std::vector<ShellID>::const_iterator
            jter = shells.begin(); jter != shells.end(); ++jter)
        {
            if(*jter == ignore1 || *jter == ignore2) continue;

            const storage_type shell(this->get_shell(*jter));
            const Real dist = boost::apply_visitor(distance_on_surf, shell);
            if(dist < radius)
            {
                for(std::size_t i=0; i<retval.size(); ++i)
                {
                    if(retval.at(i).first.first == *jter)
                    {
                        std::cerr << "Error: broken consistency in structure id\n";
                        std::cerr << "ShellContainer::list_shells_within_radius\n";
                        std::cerr << "Vertex " << pos.second << " has ...\n";
                        for(typename std::vector<VertexID>::const_iterator
                            nvi(neighborv.begin()), nve(neighborv.end());
                            nvi != nve; ++nvi)
                        {
                            std::cerr << "  VertexID = " << *nvi << std::endl;
                        }
                        std::cerr << "neighbors.\n";
                        std::cerr << "Shell " << *jter
                                  << " is found on both Vertex " << *iter
                                  << " and another Face." << std::endl;
                        face_registrator_.dump(std::cerr);
                        vertex_registrator_.dump(std::cerr);
                        assert(false);
                    }
                }
                retval.push_back(std::make_pair(
                            std::make_pair(*jter, shell), dist));
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
