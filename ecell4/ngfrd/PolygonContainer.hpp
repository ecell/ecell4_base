#ifndef ECELL4_NGFRD_POLYGON_CONTAINER_HPP
#define ECELL4_NGFRD_POLYGON_CONTAINER_HPP
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/Polygon.hpp>
#include <boost/optional.hpp>

namespace ecell4
{
namespace ngfrd
{

template<typename ObjectID> // ShellID or ParticleID
class PolygonContainer
{
public:
    explicit PolygonContainer(const std::shared_ptr<Polygon>& poly)
        : polygon_(poly)
    {}
    ~PolygonContainer() = default;
    PolygonContainer(const PolygonContainer&) = default;
    PolygonContainer(PolygonContainer&&)      = default;
    PolygonContainer& operator=(const PolygonContainer&) = default;
    PolygonContainer& operator=(PolygonContainer&&)      = default;

    boost::optional<FaceID> on_which_face(const ObjectID& id) const noexcept
    {
        const auto found = obj_to_face_.find(id);
        if(found == obj_to_face_.end())
        {
            return boost::none;
        }
        return *found;
    }
    boost::optional<EdgeID> on_which_edge(const ObjectID& id) const noexcept
    {
        const auto found = obj_to_edge_.find(id);
        if(found == obj_to_edge_.end())
        {
            return boost::none;
        }
        return *found;
    }
    boost::optional<VertexID> on_which_vertex(const ObjectID& id) const noexcept
    {
        const auto found = obj_to_vertex_.find(id);
        if(found == obj_to_vertex_.end())
        {
            return boost::none;
        }
        return *found;
    }

    boost::optional<std::vector<ObjectID> const&> objects_on(const FaceID& id) const noexcept
    {
        const auto found = face_to_obj_.find(id);
        if(found == face_to_obj_.end())
        {
            return boost::none;
        }
        return *found;
    }
    boost::optional<std::vector<ObjectID> const&> objects_on(const EdgeID& id) const noexcept
    {
        const auto found = edge_to_obj_.find(id);
        if(found == edge_to_obj_.end())
        {
            return boost::none;
        }
        return *found;
    }
    boost::optional<std::vector<ObjectID> const&> objects_on(const VertexID& id) const noexcept
    {
        const auto found = vertex_to_obj_.find(id);
        if(found == vertex_to_obj_.end())
        {
            return boost::none;
        }
        return *found;
    }

    // returns true if object is newly added.
    bool update(const ObjectID& id, const FaceID& fid)
    {
        const auto newly_added = this->remove_if_exists(id);
        if(face_to_obj_.count(fid) == 0)
        {
            face_to_obj_[fid] = std::vector<ObjectID>{};
        }
        face_to_obj_[fid].push_back(id);
        obj_to_face_[id] = fid;
        return newly_added;
    }
    bool update(const ObjectID& id, const EdgeID& eid)
    {
        const auto newly_added = this->remove_if_exists(id);
        if(edge_to_obj_.count(fid) == 0)
        {
            edge_to_obj_[fid] = std::vector<ObjectID>{};
        }
        edge_to_obj_[fid].push_back(id);
        obj_to_edge_[id] = fid;
        return newly_added;
    }
    bool update(const ObjectID& id, const VertexID& vid)
    {
        const auto newly_added = this->remove_if_exists(id);
        if(vertex_to_obj_.count(fid) == 0)
        {
            vertex_to_obj_[fid] = std::vector<ObjectID>{};
        }
        vertex_to_obj_[fid].push_back(id);
        obj_to_vertex_[id] = fid;
        return newly_added;
    }

    bool on_face  (const ObjectID& id) const noexcept {return obj_to_face_  .count(id) != 0;}
    bool on_edge  (const ObjectID& id) const noexcept {return obj_to_edge_  .count(id) != 0;}
    bool on_vertex(const ObjectID& id) const noexcept {return obj_to_vertex_.count(id) != 0;}

    void remove(const ObjectID& id)
    {
        if(const auto fid = this->on_which_face  (id)) {return remove(id, *fid);}
        if(const auto vid = this->on_which_vertex(id)) {return remove(id, *vid);}
        if(const auto eid = this->on_which_edge  (id)) {return remove(id, *eid);}
        throw_exception<NotFound>("ngfrd::PolygonContainer::remove(", id, ")");
    }

    void remove(const ObjectID& id, const FaceID&   fid)
    {
        using std::swap;
        auto& objs = face_to_obj_[fid];
        const auto found = std::find(objs.begin(), objs.end(), id);
        assert(found != objs.end());
        swap(*found, objs.back());
        objs.pop_back();

        obj_to_face_.erase(id);
        return;
    }
    void remove(const ObjectID& id, const EdgeID&   eid)
    {
        using std::swap;
        auto& objs = edge_to_obj_[fid];
        const auto found = std::find(objs.begin(), objs.end(), id);
        assert(found != objs.end());
        swap(*found, objs.back());
        objs.pop_back();

        obj_to_edge_.erase(id);
        return;
    }
    void remove(const ObjectID& id, const VertexID& vid)
    {
        using std::swap;
        auto& objs = vertex_to_obj_[fid];
        const auto found = std::find(objs.begin(), objs.end(), id);
        assert(found != objs.end());
        swap(*found, objs.back());
        objs.pop_back();

        obj_to_vertex_.erase(id);
        return;
    }

    bool diagnosis() const // for debug
    {
        bool is_ok = true;

        // --------------------------------------------------------------------
        // check obj_to_face -> face_to_obj

        std::unordered_set<ObjectID> found; // check uniqueness
        for(const auto kv : obj_to_face_)
        {
            const auto oid = kv.first;
            const auto fid = kv.second;

            if(face_to_obj_.count(fid) == 0)
            {
                std::cerr << "object " << oid << " is on " << fid << " but "
                          << " the face does not have any object" << std::endl;
                is_ok = false;
            }
            else
            {
                const auto& objs = face_to_obj_.at(fid);
                if(std::find(objs.begin(), objs.end(), oid) == objs.end())
                {
                    std::cerr << "object " << oid << " that is on " << fid
                              << " is not listed in " << fid << std::endl;
                    is_ok = false;
                }
            }

            const auto inserted = found.insert(oid);
            if(!inserted.second)
            {
                std::cerr << "object " << oid << " is assigned twice" << std::endl;
                is_ok = false;
            }
        }
        for(const auto kv : obj_to_edge_)
        {
            const auto oid = kv.first;
            const auto eid = kv.second;

            if(edge_to_obj_.count(eid) == 0)
            {
                std::cerr << "object " << oid << " is on " << eid << " but "
                          << " the edge does not have any object" << std::endl;
                is_ok = false;
            }
            else
            {
                const auto& objs = edge_to_obj_.at(eid);
                if(std::find(objs.begin(), objs.end(), oid) == objs.end())
                {
                    std::cerr << "object " << oid << " that is on " << eid
                              << " is not listed in " << eid << std::endl;
                    is_ok = false;
                }
            }

            const auto inserted = found.insert(oid);
            if(!inserted.second)
            {
                std::cerr << "object " << oid << " is assigned twice" << std::endl;
                is_ok = false;
            }
        }
        for(const auto kv : obj_to_vertex_)
        {
            const auto oid = kv.first;
            const auto vid = kv.second;

            if(vertex_to_obj_.count(vid) == 0)
            {
                std::cerr << "object " << oid << " is on " << vid << " but "
                          << " the vertex does not have any object" << std::endl;
                is_ok = false;
            }
            else
            {
                const auto& objs = vertex_to_obj_.at(vid);
                if(std::find(objs.begin(), objs.end(), oid) == objs.end())
                {
                    std::cerr << "object " << oid << " that is on " << vid
                              << " is not listed in " << vid << std::endl;
                    is_ok = false;
                }
            }
            const auto inserted = found.insert(oid);
            if(!inserted.second)
            {
                std::cerr << "object " << oid << " is assigned twice" << std::endl;
                is_ok = false;
            }
        }

        // --------------------------------------------------------------------
        // check all the objects listed in each structural element are on that
        // structure (face_to_obj -> obj_to_face)

        for(const auto& kv : face_to_obj_)
        {
            const auto fid = kv.first;
            for(const auto& oid : kv.second)
            {
                if(obj_to_face_.count(oid) == 0 || obj_to_face_.at(oid) != fid)
                {
                    std::cerr << "object " << oid << " is listed in " << fid
                              << " but not on " << fid << std::endl;
                    is_ok = false;
                }
            }
        }
        for(const auto& kv : edge_to_obj_)
        {
            const auto eid = kv.first;
            for(const auto& oid : kv.second)
            {
                if(obj_to_edge_.count(oid) == 0 || obj_to_edge_.at(oid) != eid)
                {
                    std::cerr << "object " << oid << " is listed in " << eid
                              << " but not on " << eid << std::endl;
                    is_ok = false;
                }
            }
        }
        for(const auto& kv : vertex_to_obj_)
        {
            const auto vid = kv.first;
            for(const auto& oid : kv.second)
            {
                if(obj_to_vertex_.count(oid) == 0 || obj_to_vertex_.at(oid) != vid)
                {
                    std::cerr << "object " << oid << " is listed in " << vid
                              << " but not on " << vid << std::endl;
                    is_ok = false;
                }
            }
        }
        return is_ok;
    }

private:

    bool remove_if_exists(const ObjectID& id)
    {
        if(const auto fid = on_which_face  (id)) {remove(id, *fid); return false;}
        if(const auto vid = on_which_vertex(id)) {remove(id, *vid); return false;}
        if(const auto eid = on_which_edge  (id)) {remove(id, *eid); return false;}
        return true;
    }

private:
    std::shared_ptr<Polygon> polygon_;

    std::unordered_map<ObjectID, FaceID  > obj_to_face_;
    std::unordered_map<ObjectID, EdgeID  > obj_to_edge_;
    std::unordered_map<ObjectID, VertexID> obj_to_vertex_;

    std::unordered_map<FaceID  , std::vector<ObjectID>> face_to_obj_;
    std::unordered_map<EdgeID  , std::vector<ObjectID>> edge_to_obj_;
    std::unordered_map<VertexID, std::vector<ObjectID>> vertex_to_obj_;
};

} // ngfrd
} // ecell4
#endif// ECELL4_NGFRD_POLYGON_CONTAINER_HPP
