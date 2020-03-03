#ifndef ECELL4_CORE_PERIODIC_RTREE_HPP
#define ECELL4_CORE_PERIODIC_RTREE_HPP

#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Integer3.hpp>
#include <ecell4/core/AABB.hpp>
#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/get_mapper_mf.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <utility>
#include <vector>
#include <limits>
#include <sstream>
#include <algorithm>

#ifdef WITH_HDF5
#include <ecell4/core/ParticleSpaceHDF5Writer.hpp>
#endif

namespace ecell4
{

// AABBGetter should be a functor that calculate AABB from an instance of Object
// with the following interface.
// ```
// struct AABBGetter {
//     AABB operator()(const Object& obj, const Real margin) const noexcept;
// };
// ```
// Margin is used to reduce the cost of update. Removing and inserting an object
// costs high because it requires all the nodes on top of the leaf node to be
// adjusted. It sometimes causes node re-insertion. To avoid frequent update
// while the simulation that often causes fast oscillation of the object, the
// node AABB is expanded a bit and the margin buffers those small displacements.
// Essentially the value of the margin can be set any positive value. The size
// of the margin does not affect on the accuracy but only on the runtime
// efficiency. Thus AABBGetter can be implemented in several ways, e.g.
//  1. It can completely ignore the margin.
//  2. It can add constant offset to the box size.
//  3. It can combine object information and the margin. If a object has its
//     diffusion coefficient `D`, the margin can be set as a factor of D.
// The only requirements for AABBGetter is that when margin == 0, the AABB
// should be exact (no margin).
//
// The default values of MinEntry and MaxEntry are not fine-tuned, so later
// we may need to tune it based on a benchmark result.
template<typename ObjectID, typename Object, typename AABBGetter,
         std::size_t MinEntry = 3, std::size_t MaxEntry = 8>
class PeriodicRTree
{
public:
    static_assert(!std::is_same<ObjectID, std::size_t>::value,
                  "Since node uses std::size_t as a internal node elements, "
                  "ObjectID should be distinguishable from std::size_t.");

    using box_type              = AABB;
    using box_getter_type       = AABBGetter;
    using rtree_value_type      = std::pair<box_type, std::size_t>;
    using value_type            = std::pair<ObjectID, Object>;
    using container_type        = std::vector<value_type>;
    using iterator              = typename container_type::iterator;
    using const_iterator        = typename container_type::const_iterator;
    using key_to_value_map_type = typename utils::get_mapper_mf<ObjectID, std::size_t>::type;

    static constexpr std::size_t min_entry = MinEntry;
    static constexpr std::size_t max_entry = MaxEntry;
    static constexpr std::size_t nil = std::numeric_limits<std::size_t>::max();

    struct node_type
    {
        using internal_entry_type = boost::container::static_vector<std::size_t, max_entry>;
        using leaf_entry_type     = boost::container::static_vector<ObjectID,    max_entry>;
        using entry_type          = boost::variant<internal_entry_type, leaf_entry_type>;

        struct size_invoker: boost::static_visitor<std::size_t>
        {
            template<typename T>
            std::size_t operator()(const T& v) const noexcept {return v.size();}
        };

        node_type(const internal_entry_type& inode, const std::size_t parent_,
                  const box_type& b)
            : parent(parent_), entry(inode), box(b)
        {}
        node_type(const leaf_entry_type& leaf, const std::size_t parent_,
                  const box_type& b)
            : parent(parent_), entry(leaf), box(b)
        {}
        ~node_type() = default;
        node_type(node_type const&) = default;
        node_type(node_type &&)     = default;
        node_type& operator=(node_type const&) = default;
        node_type& operator=(node_type &&)     = default;

        bool is_leaf() const noexcept {return entry.which() == 1;}
        bool empty() const noexcept
        {
            return boost::apply_visitor(size_invoker(), this->entry) == 0;
        }
        bool has_enough_storage() const noexcept
        {
            return boost::apply_visitor(size_invoker(), this->entry) < max_entry;
        }
        bool has_enough_entry()   const noexcept
        {
            return min_entry <= boost::apply_visitor(size_invoker(), this->entry);
        }

        std::size_t size() const noexcept
        {
            return boost::apply_visitor(size_invoker(), this->entry);
        }

        internal_entry_type const& inode_entry() const {return boost::get<internal_entry_type>(entry);}
        internal_entry_type&       inode_entry()       {return boost::get<internal_entry_type>(entry);}
        leaf_entry_type const&     leaf_entry()  const {return boost::get<leaf_entry_type>(entry);}
        leaf_entry_type&           leaf_entry()        {return boost::get<leaf_entry_type>(entry);}

        std::size_t parent;
        entry_type  entry;
        box_type    box;
    };

    using internal_entry_type = typename node_type::internal_entry_type;
    using leaf_entry_type     = typename node_type::leaf_entry_type;
    using tree_type           = std::vector<node_type>;
    using index_buffer_type   = std::vector<std::size_t>;

public:

    explicit PeriodicRTree(const Real3& edge_lengths, const Real margin = 0.0)
        : root_(nil), margin_(margin), edge_lengths_(edge_lengths)
    {
        // do nothing
    }
    ~PeriodicRTree() = default;
    PeriodicRTree(PeriodicRTree const&) = default;
    PeriodicRTree(PeriodicRTree &&)     = default;
    PeriodicRTree& operator=(PeriodicRTree const&) = default;
    PeriodicRTree& operator=(PeriodicRTree &&)     = default;

    std::size_t size() const noexcept
    {
        return this->container_.size();
    }

    bool empty() const noexcept {return this->root_ == nil;}
    void clear()
    {
        this->root_ = nil;
        this->tree_.clear();
        this->container_.clear();
        this->rmap_.clear();
        this->overwritable_nodes_.clear();
        return;
    }

    bool has(const ObjectID& id) const
    {
        return rmap_.count(id) != 0;
    }
    value_type get(const ObjectID& id) const
    {
        assert(id == container_.at(rmap_.at(id)).first);
        return container_.at(rmap_.at(id));
    }

    // collect values when the filter returns true.
    //
    // For example, see ParticleContainerRTreeImpl.
    // ```cpp
    // struct QueryFilter
    // {
    //     // skip based on some attribute values in a value.
    //     // If it matches, you can put an additional information to the result.
    //     boost::optional<UserDefinedInfo>
    //     operator()(const std::pair<ObjectID, Object>&, const Real3& edges);
    //
    //     // to skip non-interacting nodes purely geometric criteria.
    //     bool operator()(const AABB&,       const Real3& edges);
    // }
    // ```
    template<typename F, typename OutputIterator>
    OutputIterator query(F matches, OutputIterator out) const
    {
        if(this->empty())
        {
            return out;
        }
        return query_recursive(root_, matches, out);
    }

    // update an object. If the object corresponds to id already exists,
    // it reshapes the tree structure.
    //
    // When it reshape the structure, it considers the margin. If the node AABB
    // covers the tight (i.e. margin == 0) object AABB, we don't need to reshape
    // it. When the tight AABB of the object sticks out of the node AABB, we
    // need to re-shape it.
    bool update(const ObjectID& id, const Object& obj)
    {
        if(!this->has(id))
        {
            this->insert(id, obj);
            return true;
        }

        const auto value_idx = this->rmap_.at(id);
        const auto found = this->find_leaf(this->root_,
                                           this->container_.at(value_idx));
        if(!found)
        {
            throw std::logic_error("[error] internal error in PeriodicRTree");
        }
        const auto node_idx = found->first;
        assert(*(found->second) == id);

        const auto tight_box = this->box_getter_(obj, /*margin = */ 0.0);
        const auto& node_box = this->node_at(node_idx).box;

        if(this->is_inside(tight_box, /* is inside of */ node_box))
        {
            // the updated object is inside of the node AABB.
            // We don't need to expand the box. just replace it.
            this->container_.at(value_idx).second = obj;
        }
        else
        {
            // The updated object sticks out of the node.
            // We need to expand the box and update the tree.
            this->erase_impl(found->first, found->second);
            this->insert(id, obj);
        }
        return false;
    }
    bool update(const value_type& v)
    {
        return this->update(v.first, v.second);
    }

    void insert(const ObjectID& id, const Object& obj)
    {
        this->insert(std::make_pair(id, obj));
    }
    void insert(const value_type& v)
    {
        this->add_value(v);
        const auto box = this->box_getter_(v.second, this->margin_);
        const auto L   = this->choose_leaf(box);
        assert(node_at(L).is_leaf());

        if(node_at(L).has_enough_storage())
        {
            if(node_at(L).empty())
            {
                node_at(L).box = box;
            }
            else
            {
                node_at(L).box = this->expand(node_at(L).box, box);
            }
            this->node_at(L).leaf_entry().push_back(v.first);
            this->adjust_tree(L);
            assert(this->diagnosis());
        }
        else // the most appropreate node is already full. split it.
        {
            const auto LL = this->add_node(this->split_leaf(L, v.first, box));
            assert(L != LL);
            this->adjust_tree(L, LL);
            assert(this->diagnosis());
        }
        return ;
    }

    void erase(const ObjectID& id)
    {
        this->erase(id, this->container_.at(this->rmap_.at(id)));
        return;
    }
    void erase(const ObjectID& id, const Object& obj)
    {
        this->erase(std::make_pair(id, obj));
        return;
    }
    void erase(const value_type& v)
    {
        if(this->root_ == nil)
        {
            throw NotFound("PeriodicRTree::erase: tree is empty.");
        }

        // find_leaf returns pair{leaf_node_idx, iterator of leaf_node.entry}
        if(const auto found = this->find_leaf(this->root_, v))
        {
            this->erase_impl(found->first, found->second);
            assert(this->diagnosis());
            return;
        }

        throw NotFound("PeriodicRTree::erase: value not found.");
    }

    Real3 const& edge_lengths() const noexcept {return edge_lengths_;}
    Real3&       edge_lengths()       noexcept {return edge_lengths_;}

    // user-defined AABBGetter may be stateful (unlikely).
    AABBGetter const& get_aabb_getter() const noexcept {return box_getter_;}

    // box margin.
    Real margin() const noexcept {return margin_;}

    // container = std::vector<std::pair<ObjectID, Object>>.
    container_type const& list_objects() const {return container_;}

    // check the tree structure and relationships between nodes
    bool diagnosis() const
    {
        // -------------------------------------------------------------------
        // check number of active internal nodes
        std::size_t num_inodes  = 1; // +1 for the root
        for(std::size_t i=0; i<tree_.size(); ++i)
        {
            if(!this->is_valid_node_index(i)){continue;}

            if(!this->node_at(i).is_leaf())
            {
                num_inodes += this->node_at(i).inode_entry().size();
            }
        }
        if(tree_.size() - overwritable_nodes_.size() != num_inodes)
        {
            std::cerr << "tree_.size() is not consistent with number of "
                         "internal nodes" << std::endl;
            return false;
        }

        // -------------------------------------------------------------------
        // check the parent of a node exists and the parent contains the node
        bool root_found = false;
        for(std::size_t i=0; i<tree_.size(); ++i)
        {
            if(!this->is_valid_node_index(i)){continue;}

            if(this->node_at(i).parent == nil)
            {
                assert(!root_found);
                root_found = true;
            }
            else
            {
                const auto& e = this->node_at(this->node_at(i).parent).inode_entry();
                assert(this->is_valid_node_index(this->node_at(i).parent));
                assert(!this->node_at(this->node_at(i).parent).is_leaf());
                assert(std::find(e.begin(), e.end(), i) != e.end());
            }
        }

        // -------------------------------------------------------------------
        // check all the centroid of all the node AABBs are within the boundary
        for(std::size_t i=0; i<tree_.size(); ++i)
        {
            if(!this->is_valid_node_index(i)){continue;}

            const auto& box = this->node_at(i).box;
            const auto center = (box.upper() + box.lower()) * 0.5;
            if(!is_inside_of_boundary(center))
            {
                std::cerr << "box: " << box.lower() << ":" << box.upper() << std::endl;
                std::cerr << "center: " << center << std::endl;
            }
            assert(this->is_inside_of_boundary(center));
        }

        // -------------------------------------------------------------------
        // check all the ancester node AABBs covers the child AABBs
        for(std::size_t i=0; i<tree_.size(); ++i)
        {
            if(!this->is_valid_node_index(i)) {continue;}
            if(!node_at(i).is_leaf()) {continue;}

            if(!diagnosis_rec(i))
            {
                std::cerr << "node " << i << " is not covered by its ancestors"
                          << std::endl;
                this->dump_from_leaf_to_root(i, std::cerr);
                this->dump(std::cerr);
                return false;
            }
        }

        // -------------------------------------------------------------------
        // check consistency between id->index map and the tree
        for(std::size_t i=0; i<container_.size(); ++i)
        {
            if(rmap_.count(container_.at(i).first) != 0)
            {
                if(rmap_.at(container_.at(i).first) != i)
                {
                    std::cerr << "Object index in a container and reverse-map "
                                 "are not consistent" << std::endl;
                    return false;
                }
            }
            if(!this->find_leaf(this->root_, container_.at(i)))
            {
                std::cerr << "Object " << container_.at(i).first << ":"
                          << container_.at(i).second << " is not found in the "
                          << "tree" << std::endl;
                return false;
            }
        }
        return true;
    }

    // -------------------------------------------------------------------
    // dump parent-child relationship
    void dump(std::ostream& os) const
    {
        std::ostringstream oss;
        if(this->is_valid_node_index(this->root_)) {oss << ' ';} else {oss << '!';}
        oss << '(' << std::setw(3) << this->root_ << ") ";
        this->dump_rec(this->root_, os, oss.str());
        os << std::flush;
        return ;
    }

private:

    // -------------------------------------------------------------------
    // check all the ancester node AABBs covers the child AABBs recursively
    bool diagnosis_rec(std::size_t node_idx) const
    {
        while(node_at(node_idx).parent != nil)
        {
            const auto& node = this->node_at(node_idx);
            if(!(this->is_inside(node.box, node_at(node.parent).box)))
            {
                std::cerr << "parent  AABB: " << node_at(node.parent).box.lower()
                          << " -> " << node_at(node.parent).box.upper()
                          << std::endl;
                std::cerr << "current AABB: " << node.box.lower()
                          << " -> " << node.box.upper()
                          << std::endl;
                return false;
            }
            node_idx = node.parent;
        }
        return true;
    }

    // -------------------------------------------------------------------
    // dump parent-child relationship
    void dump_rec(const std::size_t node_idx, std::ostream& os, std::string prefix) const
    {
        const auto& node = tree_.at(node_idx);
        if(node.is_leaf())
        {
            for(const auto& entry_id : node.leaf_entry())
            {
                os << prefix;
                os << '[' << std::setw(3) << entry_id << "]\n";
            }
            return;
        }
        else
        {
            for(const std::size_t entry : node.inode_entry())
            {
                std::ostringstream oss;
                if(!this->is_valid_node_index(entry)){oss << '!';} else {oss << ' ';}
                oss << '(' << std::setw(3) << entry << ") ";
                dump_rec(entry, os, prefix + oss.str());
            }
        }
        return;
    }

    void dump_from_leaf_to_root(std::size_t node_idx, std::ostream& os) const
    {
        std::ostringstream relation;
        std::ostringstream boxes;
        do
        {
            relation << node_idx;
            const auto& node = this->node_at(node_idx);
            if(node.parent != nil) {relation << " -> ";}

            boxes << node.box.lower() << ":" << node.box.upper() << '\n';
            node_idx   = node.parent;
        }
        while(this->node_at(node_idx).parent != nil);
        os << relation.str() << std::endl;
        os << boxes.str() << std::endl;
        return;
    }

private:

    // ------------------------------------------------------------------------
    // recursively search nodes. if the node is a leaf, check values inside it.

    template<typename F, typename OutputIterator>
    OutputIterator query_recursive(const std::size_t node_idx,
            const F& matches, OutputIterator& out) const
    {
        const node_type& node = this->node_at(node_idx);
        if(node.is_leaf())
        {
            for(const auto& entry : node.leaf_entry())
            {
                const auto& value = container_.at(rmap_.at(entry));
                if(const auto info = matches(value, this->edge_lengths_))
                {
                    *out = *info;
                    ++out;
                }
            }
            return out;
        }
        // internal node. search recursively...
        for(const std::size_t entry : node.inode_entry())
        {
            const auto& node_aabb = node_at(entry).box;
            if(matches(node_aabb, this->edge_lengths_))
            {
                this->query_recursive(entry, matches, out);
            }
        }
        return out;
    }

private:

    // ------------------------------------------------------------------------
    // get node. In the debug mode (w/o -DNDEBUG), it checks the requrested node
    // is available.
    node_type&       node_at(const std::size_t i)
    {
        assert(this->is_valid_node_index(i));
        return tree_.at(i);
    }
    node_type const& node_at(const std::size_t i) const
    {
        assert(this->is_valid_node_index(i));
        return tree_.at(i);
    }

    // Note: In the original R-Tree, i.e. without periodic boundary condition,
    //       we don't need to adjust tree after we remove an object. It is
    //       because node AABB always shrink after erasing its entry. But under
    //       the periodic condition, AABB may slide to reduce its size.
    //           Let's consider the following case. When we remove the 3rd
    //       object, the node AABB slides to the left side. By removing the
    //       object, the region that was not covered by the node AABB, between
    //       object 1 and 2, will be covered.
    //
    //       .--------- boundary --------.      .--------- boundary --------.
    //       :  _____       _   ___    _ :  =>  :  _____       _          _ :
    //       : |  1  |     |2| | 3 |  |4|:  =>  : |  1  |     |2|        |4|:
    //       : |_____|     |_| |___|  |_|:  =>  : |_____|     |_|        |_|:
    //       --------]     [--node AABB---  =>  --node AABB-----]        [---
    //
    //           We need to choose a strategy to handle this problem. There can
    //       be several possible ways. First, we can keep the original AABB
    //       after removing an object. Apparently, since the older AABB covers
    //       all the child objects, we can keep the AABB. But it enlarges AABB
    //       and makes the whole R-Tree inefficient. Second, we can adjust all
    //       the ancester nodes. It might also causes inefficiency because it
    //       requires additional calculation when erasing an object. But it make
    //       node AABB smaller and the total efficiency can increase compared to
    //       the former solution.
    void erase_impl(const std::size_t node_idx,
        const typename node_type::leaf_entry_type::const_iterator value_iter)
    {
        const auto value_id = *value_iter;
        this->node_at(node_idx).leaf_entry().erase(value_iter);
        this->erase_value(rmap_.at(value_id));

        this->condense_box(this->node_at(node_idx));
        this->adjust_tree(node_idx);
        this->condense_leaf(node_idx);
        return;
    }

    // ========================================================================
    // below: construct and manage RTree structure.

    // It choose leaf object to contain the entry. It is an implementation of
    // the quadratic algorithm introduced in the paper by Guttman A. (1984)
    std::size_t choose_leaf(const box_type& entry)
    {
        // if it's empty, the entry will be inserted in a root node.
        if(this->root_ == nil)
        {
            node_type n(leaf_entry_type{/*empty*/}, /*parent = */ nil, entry);
            this->root_ = this->add_node(n);
            return this->root_;
        }

        // choose a leaf where entry will be inserted.
        //
        // search node that can contain the new entry with minimum expansion
        // until it found a leaf.
        std::size_t node_idx = this->root_;
        while(!(this->node_at(node_idx).is_leaf()))
        {
            // find the node that can cover the entry with minimum expansion
            Real diff_area_min = std::numeric_limits<Real>::max();
            Real area_min      = std::numeric_limits<Real>::max();

            // look all the child nodes and find the node that can contain the
            // new entry with minimum expansion
            const auto& node = this->node_at(node_idx);
            for(const std::size_t i : node.inode_entry())
            {
                const auto& current_box = this->node_at(i).box;

                const Real area_initial  = this->area(current_box);
                const Real area_expanded = this->area(this->expand(current_box, entry));

                const Real diff_area = area_expanded - area_initial;
                if((diff_area <  diff_area_min) ||
                   (diff_area == diff_area_min  && area_expanded < area_min))
                {
                    node_idx      = i;
                    diff_area_min = diff_area;
                    area_min      = std::min(area_min, area_expanded);
                }
            }
        }
        return node_idx;
    }

    // It expands AABBs of all the ancester nodes to make sure that
    // all the ancester nodes of `node_idx` covers the node of `node_idx`.
    void adjust_tree(std::size_t node_idx)
    {
        while(this->node_at(node_idx).parent != nil)
        {
            const auto& node   = this->node_at(node_idx);
            auto&       parent = this->node_at(node.parent);

            // if node.box is already inside of parent.box, then we don't need
            // to expand node AABBs.
            if(this->is_inside(node.box, parent.box, 0.0))
            {
                break;
            }
            parent.box = this->expand(parent.box, node.box);
            node_idx   = node.parent;
        }
        return;
    }

    // It adjusts node AABBs and add a new node to proper location.
    void adjust_tree(const std::size_t N, const std::size_t NN)
    {
        assert(N != NN);

        // we hit the root. to assign a new node, we need to make tree deeper.
        if(this->node_at(N).parent == nil)
        {
            node_type new_root(internal_entry_type{N, NN}, /*parent = */ nil,
                               this->expand(node_at(N).box, node_at(NN).box));
            this->root_  = this->add_node(std::move(new_root));

            this->node_at( N).parent = this->root_;
            this->node_at(NN).parent = this->root_;
            return;
        }
        else
        {
            const auto& node    = this->node_at(N);
            const auto& partner = this->node_at(NN);
            assert(node.parent == partner.parent);
            assert(!node_at(node.parent).is_leaf());

            const auto parent_idx = node.parent;
            auto& parent = this->node_at(parent_idx);
            parent.box = this->expand(parent.box, node.box);

            if(parent.has_enough_storage())
            {
                parent.box = this->expand(parent.box, partner.box); // assign NN
                parent.inode_entry().push_back(NN);

                // NN is assigned to this node.
                // expand AABBs of parent nodes if needed.
                return this->adjust_tree(parent_idx);
            }
            else
            {
                // NN cannot be assigned to this node.
                // split node and rearrange the tree.
                //
                // parent -+-   N }- MaxEntry
                //         +- ... }
                //         +- (NN)
                //
                //         |
                //         v
                //
                // -+-parent -+-   N }- less than MaxEntry
                //  |         +- ... }
                //  |
                //  +-PP     -+- ... }- less than MaxEntry
                //            +-  NN }
                //
                const auto PP = this->split_node(parent_idx, NN);
                assert(parent_idx != PP);
                return this->adjust_tree(parent_idx, PP);
            }
        }
    }

    // split internal nodes by the quadratic algorithm introduced by Guttman, A. (1984)
    std::size_t split_node(const std::size_t P, const std::size_t NN)
    {
        // P -+-   N }- MaxEntry
        //    +- ... }
        //    +- (NN)
        //
        //         | split into
        //         v
        //
        // -+-P  -+-   N
        //  |     +- ...
        //  |
        //  +-PP -+- ...
        //        +-  NN

        const std::size_t PP = this->add_node(
            node_type(internal_entry_type{}, node_at(P).parent, box_type()));

        assert(P  != PP);
        assert(NN != PP);

        node_type& node    = this->node_at(P);
        node_type& partner = this->node_at(PP);
        // these are internal nodes.
        assert(!node.is_leaf());
        assert(!partner.is_leaf());

        boost::container::static_vector<
            std::pair<std::size_t, box_type>, max_entry + 1> entries;
        entries.emplace_back(NN, node_at(NN).box);

        for(const auto& entry_idx : node.inode_entry())
        {
            entries.emplace_back(entry_idx, this->node_at(entry_idx).box);
        }
        node   .inode_entry().clear();
        partner.inode_entry().clear();

        // assign first 2 entries to node and partner
        {
            const auto seeds = this->pick_seeds(entries);
            assert(seeds[0] != seeds[1]);

            node   .inode_entry().push_back(entries.at(seeds[0]).first);
            partner.inode_entry().push_back(entries.at(seeds[1]).first);

            this->node_at(entries.at(seeds[0]).first).parent = P;
            this->node_at(entries.at(seeds[1]).first).parent = PP;

            node   .box = entries.at(seeds[0]).second;
            partner.box = entries.at(seeds[1]).second;

            // Remove them from entries pool. the order should be kept.
            //
            // Remove the one corresponds to the larger index first. Otherwise,
            // after removing one, index of the other one would be changed.
            entries.erase(entries.begin() + std::max(seeds[0], seeds[1]));
            entries.erase(entries.begin() + std::min(seeds[0], seeds[1]));
        }

        while(!entries.empty())
        {
            // If we need all the rest of entries to achieve min_entry,
            // use all of them.
            if(min_entry > node.size() &&
               min_entry - node.size() >= entries.size())
            {
                for(const auto idx_box : entries)
                {
                    this->node_at(idx_box.first).parent = P;
                    node.inode_entry().push_back(idx_box.first);
                    node.box = this->expand(node.box, idx_box.second);
                }
                return PP;
            }
            if(min_entry > partner.size() &&
               min_entry - partner.size() >= entries.size())
            {
                for(const auto idx_box : entries)
                {
                    this->node_at(idx_box.first).parent = PP;
                    partner.inode_entry().push_back(idx_box.first);
                    partner.box = this->expand(partner.box, idx_box.second);
                }
                return PP;
            }

            // choose which entry will be assigned to which node
            const auto next = this->pick_next(entries, node.box, partner.box);
            if(next.second)
            {
                node.inode_entry().push_back(entries.at(next.first).first);
                this->node_at(entries.at(next.first).first).parent = P;
                node.box = this->expand(node.box, entries.at(next.first).second);
            }
            else
            {
                partner.inode_entry().push_back(entries.at(next.first).first);
                this->node_at(entries.at(next.first).first).parent = PP;
                partner.box = this->expand(partner.box, entries.at(next.first).second);
            }
            entries.erase(entries.begin() + next.first);
        }
        this->node_at(P)  = node;
        this->node_at(PP) = partner;
        return PP;
    }

    // split leaf nodes by the quadratic algorithm introduced by Guttman, A. (1984)
    node_type split_leaf(const std::size_t N, const ObjectID& vid,
                         const box_type& entry)
    {
        node_type& node = node_at(N);
        node_type  partner(leaf_entry_type{}, node.parent, box_type());
        assert(node.is_leaf());

        boost::container::static_vector<
            std::pair<ObjectID, box_type>, max_entry + 1> entries;
        entries.push_back(std::make_pair(vid, entry));

        for(const auto& entry_id : node.leaf_entry())
        {
            entries.emplace_back(entry_id,
                box_getter_(container_.at(rmap_.at(entry_id)).second, margin_));
        }

        node   .leaf_entry().clear();
        partner.leaf_entry().clear();

        // assign first 2 entries to node and partner
        {
            const auto seeds = this->pick_seeds(entries);
            node   .leaf_entry().push_back(entries.at(seeds[0]).first);
            partner.leaf_entry().push_back(entries.at(seeds[1]).first);

            node   .box = entries.at(seeds[0]).second;
            partner.box = entries.at(seeds[1]).second;

            // remove them from entries pool
            //
            // Remove the one corresponds to the larger index first. Otherwise,
            // after removing one, index of the other one would be changed.
            entries.erase(entries.begin() + std::max(seeds[0], seeds[1]));
            entries.erase(entries.begin() + std::min(seeds[0], seeds[1]));
        }

        while(!entries.empty())
        {
            if(min_entry > node.size() &&
               min_entry - node.size() >= entries.size())
            {
                for(const auto& idx_box : entries)
                {
                    node.leaf_entry().push_back(idx_box.first);
                    node.box = this->expand(node.box, idx_box.second);
                }
                return partner;
            }
            if(min_entry > partner.size() &&
               min_entry - partner.size() >= entries.size())
            {
                for(const auto& idx_box : entries)
                {
                    partner.leaf_entry().push_back(idx_box.first);
                    partner.box = this->expand(partner.box, idx_box.second);
                }
                return partner;
            }

            const auto next = this->pick_next(entries, node.box, partner.box);
            if(next.second) // next is for node
            {
                node.leaf_entry().push_back(entries.at(next.first).first);
                node.box = this->expand(node.box, entries.at(next.first).second);
            }
            else // next is for partner
            {
                partner.leaf_entry().push_back(entries.at(next.first).first);
                partner.box = this->expand(partner.box, entries.at(next.first).second);
            }
            entries.erase(entries.begin() + next.first);
        }
        return partner;
    }

    // auxiliary function for the quadratic algorithm.
    template<typename T>
    std::array<std::size_t, 2>
    pick_seeds(const boost::container::static_vector<
            std::pair<T, box_type>, max_entry+1>& entries) const
    {
        assert(entries.size() >= 2);

        std::array<std::size_t, 2> retval{{
            std::numeric_limits<std::size_t>::max(),
            std::numeric_limits<std::size_t>::max()
        }};

        boost::container::static_vector<Real, max_entry+1> areas;
        for(const auto& idx_box : entries)
        {
            areas.push_back(this->area(idx_box.second));
        }

        // Choose a pair that are most distant to each other.
        // When we split a node, we always start such a pair.
        Real max_d = -std::numeric_limits<Real>::max();
        for(std::size_t i=0; i+1<entries.size(); ++i)
        {
            const auto& ibox = entries.at(i).second;
            const Real iarea = areas.at(i);
            for(std::size_t j=i+1; j<entries.size(); ++j)
            {
                const auto& jbox   = entries.at(j).second;
                const auto  merged = this->expand(ibox, jbox);
                const Real  d      = this->area(merged) - iarea - areas.at(j);
                if(max_d < d)
                {
                    max_d     = d;
                    retval[0] = i;
                    retval[1] = j;
                }
            }
        }
        return retval;
    }

    // auxiliary function for the quadratic algorithm.
    template<typename T>
    std::pair<std::size_t, bool>
    pick_next(const boost::container::static_vector<
            std::pair<T, box_type>, max_entry+1>& entries,
            const box_type& node, const box_type& ptnr) const
    {
        assert(!entries.empty());

        bool is_node;
        Real max_dd = -1;
        std::size_t idx;
        for(std::size_t i=0; i<entries.size(); ++i)
        {
            const auto& idx_box = entries.at(i);
            const auto box1 = this->expand(node, idx_box.second);
            const auto box2 = this->expand(ptnr, idx_box.second);

            const Real d1 = this->area(box1) - this->area(node);
            const Real d2 = this->area(box2) - this->area(ptnr);
            const Real dd = d1 - d2;
            if(max_dd < std::abs(dd))
            {
                max_dd  = std::abs(dd);
                idx     = i;
                is_node = (dd < 0);
            }
        }
        return std::make_pair(idx, is_node);
    }


    // find leaf node that contains the entry value in a recursive manner.
    //
    // returns pairof{leaf-node-index, entry-index-in-node}
    boost::optional<std::pair<std::size_t, typename
                              node_type::leaf_entry_type::const_iterator>>
    find_leaf(std::size_t node_idx, const value_type& entry) const
    {
        const node_type& node = this->node_at(node_idx);
        const auto tight_box  = this->box_getter_(entry.second, 0.0);

        if(!(this->is_inside(tight_box, node.box)))
        {
            return boost::none;
        }
        if(node.is_leaf())
        {
            for(auto i=node.leaf_entry().begin(), e=node.leaf_entry().end(); i!=e; ++i)
            {
                // if the ID is the same, then the objects should be the same.
                if(container_.at(rmap_.at(*i)).first == entry.first)
                {
                    return std::make_pair(node_idx, i);
                }
            }
            return boost::none;
        }
        else // node is an internal node
        {
            for(auto i=node.inode_entry().begin(), e=node.inode_entry().end(); i!=e; ++i)
            {
                if(!(this->is_inside(tight_box, this->node_at(*i).box)))
                {
                    continue;
                }
                if(const auto found = this->find_leaf(*i, entry))
                {
                    return found;
                }
            }
            return boost::none;
        }
    }

    // check the number of entries in the Nth leaf node. If it has too few
    // entries, it removes the node and balance the tree.
    void condense_leaf(const std::size_t N)
    {
        const node_type& node = this->node_at(N);
        assert(node.is_leaf());

        if(node.has_enough_entry() || node.parent == nil)
        {
            return; // nothing is needed.
        }

        // Leaf has too few entries. The tree structure should be balanced.
        // Here, the naivest approach is chosen. First remove the leaf node
        // that has too few children and re-insert all of them later.

        // copy objects in the leaf node that is being removed to re-insert them later
        std::vector<value_type> eliminated_objs;
        eliminated_objs.reserve(node.size());
        for(const auto& id: node.leaf_entry())
        {
            const std::size_t idx = rmap_.at(id);
            eliminated_objs.push_back(this->container_.at(idx));
            this->erase_value(idx);
        }

        const auto parent_idx = node.parent;

        // erase the node N from its parent and condense aabb
        auto found = std::find(this->node_at(node.parent).inode_entry().begin(),
                               this->node_at(node.parent).inode_entry().end(), N);
        assert(found != this->node_at(node.parent).inode_entry().end());

        this->erase_node(*found);
        this->node_at(parent_idx).inode_entry().erase(found);

        // condense node parent box without the leaf node N.
        this->condense_box(this->node_at(parent_idx));
        this->adjust_tree(parent_idx);
        this->condense_node(parent_idx);

        // re-insert entries that were in node N
        for(const auto& obj : eliminated_objs)
        {
            this->insert(obj);
        }
        return;
    }

    // check the number of entries in the Nth internal node. If it has too few
    // entries, it removes the node and balance the tree.
    void condense_node(const std::size_t N)
    {
        const node_type& node = this->node_at(N);
        assert(!node.is_leaf());

        if(node.has_enough_entry())
        {
            return; // if the node has enough number of entries, then it's okay.
        }
        if(node.parent == nil && node.size() == 1)
        {
            // if the root has only one entry, then the child node of the current
            // root node should be the root node, no?
            this->root_ = node.inode_entry().front();
            this->erase_node(N);
            assert(!this->is_valid_node_index(N));
            return;
        }

        // collect index of nodes that are children of the node to be removed
        const std::vector<std::size_t> eliminated_nodes(
                node.inode_entry().begin(), node.inode_entry().end());
        const auto parent_idx = node.parent;

        // erase the node N from its parent and condense its aabb
        auto found = std::find(this->node_at(parent_idx).inode_entry().begin(),
                               this->node_at(parent_idx).inode_entry().end(), N);
        assert(found != this->node_at(parent_idx).inode_entry().end());

        // remove the node from its parent.entry
        this->erase_node(*found);
        this->node_at(parent_idx).inode_entry().erase(found);
        this->condense_box(this->node_at(parent_idx));

        // re-insert nodes eliminated from node N
        for(const auto& idx : eliminated_nodes)
        {
            this->re_insert(idx);
        }
        this->condense_node(parent_idx);
        return;
    }

    // re-calculate the AABB of a node to make sure that the node covers all the
    // child nodes.
    void condense_box(node_type& node) const
    {
        if(node.empty())
        {
            // If the entry is empty, the node will soon be eliminated from its
            // parent node via `condense_leaf` or `condense_node`.
            // We can delete this operation and just return from this method,
            // but I'm not completely sure yet.
            const auto center = (node.box.upper() + node.box.lower()) * 0.5;
            node.box = box_type(center, center);
            return;
        }
        if(node.is_leaf())
        {
            const auto& entries = node.leaf_entry();
            node.box = box_getter_(
                this->container_.at(rmap_.at(entries.front())).second, margin_);
            for(auto i = std::next(entries.begin()); i != entries.end(); ++i)
            {
                node.box = this->expand(node.box,
                    box_getter_(container_.at(rmap_.at(*i)).second, margin_));
            }
        }
        else
        {
            const auto& entries = node.inode_entry();
            node.box = this->node_at(entries.front()).box;
            for(auto i = std::next(entries.begin()); i != entries.end(); ++i)
            {
                node.box = this->expand(node.box, this->node_at(*i).box);
            }
        }
        return;
    }

    // re-insert nodes that are temporary removed from its (previous) parent
    // node to balance the number of nodes. This will only be called from
    // `condense_node`.
    void re_insert(const std::size_t N)
    {
        // reset connection to the parent!

        // insert node to its proper parent. to find the parent of this node N,
        // add 1 to level. root node should NOT come here.
        const std::size_t level = this->level_of(N) + 1;
        const box_type&   entry = this->node_at(N).box;
        const std::size_t L = this->choose_node_with_level(entry, level);
        assert(!this->node_at(L).is_leaf());

        if(node_at(L).has_enough_storage())
        {
            node_at(L).inode_entry().push_back(N);
            node_at(N).parent = L;
            node_at(L).box    = this->expand(this->node_at(L).box, entry);
            this->adjust_tree(L);
        }
        else
        {
            const std::size_t LL = this->split_node(L, N);
            this->adjust_tree(L, LL);
        }
        return;
    }

    // choose where the node should be inserted by using the box size and
    // the level of the node.
    std::size_t
    choose_node_with_level(const box_type& entry, const std::size_t level)
    {
        std::size_t node_idx = this->root_;
        if(this->level_of(node_idx) < level)
        {
            throw std::logic_error("PeriodicRTree: No such a high level node.");
        }

        while(this->level_of(node_idx) != level)
        {
            Real diff_area_min = std::numeric_limits<Real>::max();
            Real area_min      = std::numeric_limits<Real>::max();

            const node_type& node = this->node_at(node_idx);
            for(const auto& entry_idx : node.inode_entry())
            {
                const auto& entry_box = node_at(entry_idx).box;
                const Real area_initial = this->area(entry_box);
                const box_type      box = this->expand(entry_box, entry);

                const Real area_expanded = this->area(box);
                const Real diff_area     = area_expanded - area_initial;
                if(diff_area <  diff_area_min ||
                  (diff_area == diff_area_min && area_expanded < area_min))
                {
                    node_idx      = entry_idx;
                    diff_area_min = diff_area;
                    area_min      = std::min(area_expanded, area_min);
                }
            }
        }
        return node_idx;
    }

    // the naivest approach; going down until it hits leaf.
    // If the node is a leaf, level == 0.
    std::size_t level_of(std::size_t node_idx) const
    {
        std::size_t level = 0;
        while(!(node_at(node_idx).is_leaf()))
        {
            ++level;
            node_idx = node_at(node_idx).inode_entry().front();
        }
        return level;
    }

private:

    // ========================================================================
    // utility member methods to handle `container_` and `tree_`.

    // ------------------------------------------------------------------------
    // for nodes
    //
    // RTree stores nodes by their index in a container, `tree_`.
    // Thus the indices should be kept when we insert a new node or remove an
    // old node.

    // re-use overwritable region if it exists.
    std::size_t add_node(const node_type& n)
    {
        if(overwritable_nodes_.empty())
        {
            const std::size_t new_index = tree_.size();
            tree_.push_back(n);
            return new_index;
        }
        const std::size_t new_index = overwritable_nodes_.back();
        overwritable_nodes_.pop_back();
        node_at(new_index) = n;
        return new_index;
    }
    // mark index `i` overritable and fill old value by the default value.
    void erase_node(const std::size_t i)
    {
        node_at(i) = node_type(internal_entry_type{}, nil, box_type());
        overwritable_nodes_.push_back(i);
        assert(!this->is_valid_node_index(i));
        return;
    }
    bool is_valid_node_index(const std::size_t i) const noexcept
    {
        return std::find(overwritable_nodes_.begin(),
                         overwritable_nodes_.end(), i) ==
            overwritable_nodes_.end();
    }

    // ------------------------------------------------------------------------
    // for values

    std::size_t add_value(const value_type& v)
    {
        const std::size_t new_index = container_.size();
        container_.push_back(v);
        rmap_[v.first] = new_index;
        return new_index;
    }
    void erase_value(const std::size_t i)
    {
        using std::swap;
        swap(container_.at(i), container_.back());

        assert(rmap_.at(container_.at(i).first) == container_.size() - 1);
        rmap_[container_.at(i).first] = i;

        const auto status = rmap_.erase(container_.back().first);
        assert(status == 1);
        (void)status;

        container_.pop_back();
        return ;
    }

private:

    // ------------------------------------------------------------------------
    // geometric methods to modify AABB under periodic boundary condition.

    // upper/lower of AABB are not adjusted to the PBC but the center is aligned.
    // so the width is always correct.
    Real area(const box_type& box) const noexcept
    {
        const auto dx = box.upper() - box.lower();
        return dx[0] * dx[1] * dx[2];
    }

    // merge two AABBs under the PBC.
    //
    // It is guaranteed that the resulting AABB contains both lhs and rhs.
    box_type expand(const box_type& lhs, const box_type& rhs) const noexcept
    {
        // assuming that upper/lower can stick out of the boundary
        const auto lc = (lhs.upper() + lhs.lower()) * 0.5; // center of lhs
        const auto rc = (rhs.upper() + rhs.lower()) * 0.5; // center of rhs
        const auto rr =  rhs.upper() - rc;                 // radius of rhs
        const auto dc = restrict_direction(rc - lc); // distance between centers
        const auto l1 = lhs.lower();
        const auto u1 = lhs.upper();
        const auto l2 = lc + dc - rr; // boundary-adjusted rhs' lower bound
        const auto u2 = lc + dc + rr; // ditto

        Real3 up, lw;
        for(std::size_t i=0; i<3; ++i)
        {
            lw[i] = std::min(l1[i], l2[i]);
            up[i] = std::max(u1[i], u2[i]);
        }
        const auto c = (up + lw) * 0.5;
        const auto d = restrict_position(c) - c;

        box_type expanded(lw + d, up + d);

        assert(is_inside(lhs, expanded, 1e-8));
        assert(is_inside(rhs, expanded, 1e-8));

        return expanded;
    }

    // check if two AABBs intersects each other, under the PBC.
    bool intersects(const box_type& lhs, const box_type& rhs,
                    const Real tol = 1e-8) const noexcept
    {
        const auto lc = (lhs.upper() + lhs.lower()) * 0.5; // center of lhs
        const auto lr = (lhs.upper() - lhs.lower()) * 0.5; // radius of lhs
        const auto rc = (rhs.upper() + rhs.lower()) * 0.5; // center of rhs
        const auto rr = (rhs.upper() - rhs.lower()) * 0.5; // radius of rhs

        const auto r2 = lr + rr;
        const auto dc = ecell4::abs(this->restrict_direction(lc - rc));

        // if they are close along all the axes, they intersects each other.
        return ((dc[0] - r2[0]) <= tol) &&
               ((dc[1] - r2[1]) <= tol) &&
               ((dc[2] - r2[2]) <= tol);
    }

    // if lhs is inside of rhs, return true.
    bool is_inside(const box_type& lhs, const box_type& rhs,
                   const Real tol = 1e-8) const noexcept
    {
        const auto lc = (lhs.upper() + lhs.lower()) * 0.5; // center of lhs
        const auto lr = (lhs.upper() - lhs.lower()) * 0.5; // radius of lhs
        const auto rc = (rhs.upper() + rhs.lower()) * 0.5; // center of rhs
        const auto rr = (rhs.upper() - rhs.lower()) * 0.5; // radius of rhs

        const auto dr = rr - lr; // radius difference
        const auto dc = abs(this->restrict_direction(lc - rc)); // distance between centers

        // if the radius (of the right hand side) is larger than the half width,
        // that means that the right hand side wraps the whole world.
        return ((dc[0] - dr[0]) <= tol || (edge_lengths_[0] * 0.5 <= rr[0])) &&
               ((dc[1] - dr[1]) <= tol || (edge_lengths_[1] * 0.5 <= rr[1])) &&
               ((dc[2] - dr[2]) <= tol || (edge_lengths_[2] * 0.5 <= rr[2]));
    }

    bool is_inside_of_boundary(const Real3 r, const Real tol = 1e-8) const noexcept
    {
        const auto rel_tol = edge_lengths_ * tol;
        if(r[0] < -rel_tol[0] || edge_lengths_[0] + rel_tol[0] < r[0]) {return false;}
        if(r[1] < -rel_tol[1] || edge_lengths_[1] + rel_tol[1] < r[1]) {return false;}
        if(r[2] < -rel_tol[2] || edge_lengths_[2] + rel_tol[2] < r[2]) {return false;}
        return true;
    }

    Real3 restrict_position(Real3 r) const noexcept
    {
        return modulo(r, edge_lengths_);
    }
    Real3 restrict_direction(Real3 dr) const noexcept
    {
        // Assuming that ...
        //  - dr = r1 - r2
        //  - Both r1 and r2 are inside of the boundary.
        const auto halfw = edge_lengths_ * 0.5;
        if     (   dr[0] < -halfw[0]) {dr[0] += edge_lengths_[0];}
        else if(halfw[0] <=    dr[0]) {dr[0] -= edge_lengths_[0];}
        if     (   dr[1] < -halfw[1]) {dr[1] += edge_lengths_[1];}
        else if(halfw[1] <=    dr[1]) {dr[1] -= edge_lengths_[1];}
        if     (   dr[2] < -halfw[2]) {dr[2] += edge_lengths_[2];}
        else if(halfw[2] <=    dr[2]) {dr[2] -= edge_lengths_[2];}
        return dr;
    }

private:

    std::size_t           root_;                // index of the root node.
    Real                  margin_;              // margin of the Box.
    Real3                 edge_lengths_;        // (Lx, Ly, Lz) of PBC.
    tree_type             tree_;                // vector of nodes.
    container_type        container_;           // vector of Objects.
    key_to_value_map_type rmap_;                // map from ObjectID to index
    index_buffer_type     overwritable_nodes_;  // list "already-removed" nodes
    box_getter_type       box_getter_;          // AABBGetter can be stateful.
};

template<typename ObjID, typename Obj, typename Box, std::size_t MinE, std::size_t MaxE>
constexpr std::size_t PeriodicRTree<ObjID, Obj, Box, MinE, MaxE>::nil;
template<typename ObjID, typename Obj, typename Box, std::size_t MinE, std::size_t MaxE>
constexpr std::size_t PeriodicRTree<ObjID, Obj, Box, MinE, MaxE>::min_entry;
template<typename ObjID, typename Obj, typename Box, std::size_t MinE, std::size_t MaxE>
constexpr std::size_t PeriodicRTree<ObjID, Obj, Box, MinE, MaxE>::max_entry;

} // ecell4
#endif// ECELL4_PERIODIC_RTREE_IMPL_HPP
