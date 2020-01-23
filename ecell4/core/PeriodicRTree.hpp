#ifndef ECELL4_CORE_PERIODIC_RTREE_HPP
#define ECELL4_CORE_PERIODIC_RTREE_HPP

#include <ecell4/core/Real3.hpp>
#include <ecell4/core/Integer3.hpp>
#include <ecell4/core/AABB.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/optional.hpp>
#include <utility>
#include <unordered_map>
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
    using box_type              = AABB;
    using box_getter_type       = AABBGetter;
    using rtree_value_type      = std::pair<box_type, std::size_t>;
    using value_type            = std::pair<ObjectID, Object>;
    using container_type        = std::vector<value_type>;
    using key_to_value_map_type = std::unordered_map<ObjectID, std::size_t>;

    static constexpr std::size_t min_entry = MinEntry;
    static constexpr std::size_t max_entry = MaxEntry;
    static constexpr std::size_t nil = std::numeric_limits<std::size_t>::max();

    struct node_type
    {
        using entry_type =
            boost::container::static_vector<std::size_t, max_entry>;
        using iterator       = typename entry_type::iterator;
        using const_iterator = typename entry_type::const_iterator;

        node_type(const bool is_leaf_, const std::size_t parent_)
            : is_leaf(is_leaf_), parent(parent_)
        {}
        ~node_type() = default;
        node_type(node_type const&) = default;
        node_type(node_type &&)     = default;
        node_type& operator=(node_type const&) = default;
        node_type& operator=(node_type &&)     = default;

        bool has_enough_storage() const noexcept {return entry.size() < max_entry;}
        bool has_enough_entry()   const noexcept {return min_entry <= entry.size();}

        bool        is_leaf;
        std::size_t parent;
        entry_type  entry;
        box_type    box;
    };

    typedef std::vector<node_type>   tree_type;
    typedef std::vector<std::size_t> index_buffer_type;

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

    std::size_t size() const noexcept {return container_.size();}
    bool       empty() const noexcept {return this->root_ == nil;}
    void       clear()
    {
        this->root_ = nil;
        this->margin_ = 0.0;
        this->tree_.clear();
        this->container_.clear();
        this->rmap_.clear();
        this->overwritable_values_.clear();
        this->overwritable_nodes_.clear();
        return;
    }

    bool has(const ObjectID& id) const
    {
        return rmap_.count(id) != 0;
    }
    value_type get(const ObjectID& id) const
    {
        return container_.at(rmap_.at(id));
    }

    // collect values when the filter returns true.
    // ```cpp
    // struct QueryFilter {
    //     // skip based on some attribute values in a value.
    //     // return true if the value should be collected.
    //     bool operator(const std::pair<ObjectID, Object>&, const Real3& edges);
    //
    //     // to skip non-interacting nodes purely geometric criteria.
    //     bool operator(const AABB&,       const Real3& edges);
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
        const auto found = this->find_leaf(this->root_,
                this->container_.at(this->rmap_.at(id)));
        if(!found)
        {
            throw std::logic_error("[error] internal error in PeriodicRTree");
        }
        const auto node_idx  =   found->first;
        const auto value_idx = *(found->second);

        const auto tight_box = this->box_getter_(obj, 0.0);
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
            this->erase_impl(node_idx, found->second);
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
        const auto idx = this->add_value(v);
        const auto box = this->box_getter_(v.second, this->margin_);
        const auto L   = this->choose_leaf(box);

        if(node_at(L).has_enough_storage())
        {
            node_at(L).entry.push_back(idx);
            node_at(L).box = this->expand(node_at(L).box, box);
            this->adjust_tree(L);
        }
        else // the most appropreate node is already full. split it.
        {
            const auto LL = this->add_node(this->split_leaf(L, idx, box));
            assert(L != LL);
            this->adjust_tree(L, LL);
        }
        return ;
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
            throw std::out_of_range("PeriodicRTree::erase: tree is empty.");
        }

        // find_leaf returns pair{leaf_node_idx, iterator of leaf_node.entry}
        if(const auto found = this->find_leaf(this->root_, v))
        {
            this->erase_impl(found->first, found->second);
            return;
        }
        else
        {
            throw std::out_of_range("PeriodicRTree::erase: value not found.");
        }
    }

    Real3 const& edge_lengths() const noexcept {return edge_lengths_;}

    // user-defined AABBGetter may be stateful (unlikely).
    AABBGetter const& get_aabb_getter() const noexcept {return box_getter_;}

    // box margin.
    Real margin() const noexcept {return margin_;}

    // it does not return a reference, but returns a copy.
    // Since this class handles values via index, some values can be in a
    // invalid state. Before returning it, we need to remove those
    // "already-removed" values.
    std::vector<std::pair<ObjectID, Object>> list_objects() const
    {
        std::vector<std::pair<ObjectID, Object>> retval = container_;

        std::vector<std::size_t> overwritable(overwritable_values_.begin(),
                                              overwritable_values_.end());

        std::sort(overwritable.begin(), overwritable.end());
        for(auto ri=overwritable.crbegin();  ri!=overwritable.crend(); ++ri)
        {
            retval.erase(retval.begin() + *ri);
        }
        return retval;
    }

    // check the tree structure and relationships between nodes
    bool diagnosis() const
    {
        // -------------------------------------------------------------------
        // check number of active objects and internal nodes
        std::size_t num_objects = 0;
        std::size_t num_inodes  = 1; // +1 for the root
        for(std::size_t i=0; i<tree_.size(); ++i)
        {
            if(!this->is_valid_node_index(i)){continue;}

            if(this->node_at(i).is_leaf)
            {
                num_objects += this->node_at(i).entry.size();
            }
            else
            {
                num_inodes += this->node_at(i).entry.size();
            }
        }
        assert(list_objects().size() == num_objects);
        assert(tree_.size() - overwritable_nodes_.size() == num_inodes);

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
                const auto& e = this->node_at(this->node_at(i).parent).entry;
                assert(this->is_valid_node_index(this->node_at(i).parent));
                assert(!this->node_at(this->node_at(i).parent).is_leaf);
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
            assert(this->is_inside_of_boundary(center));
        }

        // -------------------------------------------------------------------
        // check all the ancester node AABBs covers the child AABBs
        for(std::size_t i=0; i<tree_.size(); ++i)
        {
            if(!this->is_valid_node_index(i)) {continue;}
            if(node_at(i).is_leaf && !diagnosis_rec(i))
            {
                return false;
            }
        }

        // -------------------------------------------------------------------
        // check consistency between id->index map and the tree
        for(std::size_t i=0; i<container_.size(); ++i)
        {
            if(!this->is_valid_value_index(i)) {continue;}

            if(rmap_.count(container_.at(i).first) != 0)
            {
                if(rmap_.at(container_.at(i).first) != i)
                {
                    return false;
                }
            }

            if(!this->find_leaf(this->root_, container_.at(i)))
            {
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
        if(node.is_leaf)
        {
            for(const std::size_t entry : node.entry)
            {
                os << prefix;
                if(!this->is_valid_value_index(entry)){os << '!';} else {os << ' ';}
                os << '[' << std::setw(3) << entry << "]\n";
            }
            return;
        }

        for(const std::size_t entry : node.entry)
        {
            std::ostringstream oss;
            if(!this->is_valid_node_index(entry)){oss << '!';} else {oss << ' ';}
            oss << '(' << std::setw(3) << entry << ") ";
            dump_rec(entry, os, prefix + oss.str());
        }
        return;
    }

private:

    // ------------------------------------------------------------------------
    // recursively search nodes. if the node is a leaf, check values inside it.

    template<typename F, typename OutputIterator>
    OutputIterator query_recursive(
            const std::size_t node_idx, F matches, OutputIterator& out) const
    {
        const node_type& node = this->node_at(node_idx);
        if(node.is_leaf)
        {
            for(const std::size_t entry : node.entry)
            {
                const auto& value = container_.at(entry);
                if(matches(value, this->edge_lengths_))
                {
                    *out = value;
                    ++out;
                }
            }
            return out;
        }
        // internal node. search recursively...
        for(const std::size_t entry : node.entry)
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

    // Note: In the original R-Tree, without periodic boundary condition,
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
                    const typename node_type::const_iterator value_idx)
    {
        this->node_at(node_idx).entry.erase(found->second);
        this->erase_value(value_idx);

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
            node_type n(/*is_leaf = */true, /*parent = */ nil);
            n.box = entry;
            this->root_ = this->add_node(n);
            return this->root_;
        }

        // choose a leaf where entry will be inserted.
        std::size_t node_idx = this->root_;
        while(!(this->node_at(node_idx).is_leaf))
        {
            // find the node that can cover the entry with minimum expansion
            Real diff_area_min = std::numeric_limits<Real>::max();
            Real area_min      = std::numeric_limits<Real>::max();

            const auto& node = this->node_at(node_idx);
            for(const std::size_t i : node.entry)
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

    // It adjusts AABBs of all the ancester nodes of `node_idx` to make sure
    // that all the ancester nodes covers the node of `node_idx`.
    void adjust_tree(std::size_t node_idx)
    {
        while(node_at(node_idx).parent != nil)
        {
            const auto& node   = this->node_at(node_idx);
            auto&       parent = this->node_at(node.parent);

            // if node.box is already inside of parent.box, then we don't need
            // to expand node AABBs.
            if(this->is_inside(node.box, parent.box))
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
        if(node_at(N).parent == nil)
        {
            node_type new_root(/*is_leaf = */false, /*parent = */ nil);
            new_root.entry.push_back(N);
            new_root.entry.push_back(NN);
            new_root.box = this->expand(node_at(N).box, node_at(NN).box);
            this->root_  = this->add_node(std::move(new_root));

            this->node_at( N).parent = this->root_;
            this->node_at(NN).parent = this->root_;
            return;
        }
        else
        {
            const auto& node    = node_at(N);
            const auto& partner = node_at(NN);
            assert(node.parent == partner.parent);

            const auto parent_idx = node.parent;
            auto& parent = node_at(parent_idx);
            parent.box = this->expand(parent.box, node.box);

            if(parent.has_enough_storage())
            {
                parent.box = this->expand(parent.box, partner.box); // assign NN
                parent.entry.push_back(NN);

                // NN is assigned to this node. expand AABBs of parent nodes
                // if needed.
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
                // -+-parent -+-   N
                //  |         +- ...
                //  +-PP     -+- ...
                //            +-  NN
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
        //         |
        //         v
        //
        // -+-P  -+-   N
        //  |     +- ...
        //  +-PP -+- ...
        //        +-  NN

        const std::size_t PP = this->add_node(node_type(false, node_at(P).parent));
        assert(P  != PP);
        assert(NN != PP);
        node_type& node    = node_at(P);
        node_type& partner = node_at(PP);

        assert(!node.is_leaf);
        assert(!partner.is_leaf);

        boost::container::static_vector<
            std::pair<std::size_t, box_type>, max_entry + 1> entries;
        entries.emplace_back(NN, node_at(NN).box);

        for(const auto& entry_idx : node.entry)
        {
            entries.emplace_back(entry_idx, node_at(entry_idx).box);
        }
        node.entry.clear();
        partner.entry.clear();

        // assign first 2 entries to node and partner
        {
            const auto seeds = this->pick_seeds(entries);
            assert(seeds[0] != seeds[1]);

            node   .entry.push_back(entries.at(seeds[0]).first);
            partner.entry.push_back(entries.at(seeds[1]).first);

            node_at(entries.at(seeds[0]).first).parent = P;
            node_at(entries.at(seeds[1]).first).parent = PP;

            node   .box = entries.at(seeds[0]).second;
            partner.box = entries.at(seeds[1]).second;

            // remove them from entries pool. the order should be kept.
            entries.erase(entries.begin() + std::max(seeds[0], seeds[1]));
            entries.erase(entries.begin() + std::min(seeds[0], seeds[1]));
        }

        while(!entries.empty())
        {
            // If we need all the rest of entries to achieve min_entry,
            // use all of them.
            if(min_entry > node.entry.size() &&
               min_entry - node.entry.size() >= entries.size())
            {
                for(const auto idx_box : entries)
                {
                    node_at(idx_box.first).parent = P;
                    node.entry.push_back(idx_box.first);
                    node.box = this->expand(node.box, idx_box.second);
                }
                return PP;
            }
            if(min_entry > partner.entry.size() &&
               min_entry - partner.entry.size() >= entries.size())
            {
                for(const auto idx_box : entries)
                {
                    node_at(idx_box.first).parent = PP;
                    partner.entry.push_back(idx_box.first);
                    partner.box = this->expand(partner.box, idx_box.second);
                }
                return PP;
            }

            // choose which entry will be assigned to which node
            const auto next = this->pick_next(entries, node.box, partner.box);
            if(next.second)
            {
                node.entry.push_back(entries.at(next.first).first);
                node_at(entries.at(next.first).first).parent = P;
                node.box = this->expand(node.box, entries.at(next.first).second);
            }
            else
            {
                partner.entry.push_back(entries.at(next.first).first);
                node_at(entries.at(next.first).first).parent = PP;
                partner.box = this->expand(partner.box, entries.at(next.first).second);
            }
            entries.erase(entries.begin() + next.first);
        }
        node_at(P)  = node;
        node_at(PP) = partner;
        return PP;
    }

    // split leaf nodes by the quadratic algorithm introduced by Guttman, A. (1984)
    node_type split_leaf(const std::size_t N, const std::size_t vidx,
                         const box_type& entry)
    {
        node_type& node = node_at(N);
        node_type  partner(true, node.parent);
        assert(node.is_leaf);

        boost::container::static_vector<
            std::pair<std::size_t, box_type>, max_entry + 1> entries;
        entries.push_back(std::make_pair(vidx, entry));

        for(const auto& entry_idx : node.entry)
        {
            entries.emplace_back(entry_idx,
                    box_getter_(container_.at(entry_idx).second, this->margin_));
        }

        node   .entry.clear();
        partner.entry.clear();

        // assign first 2 entries to node and partner
        {
            const auto seeds = this->pick_seeds(entries);
            node   .entry.push_back(entries.at(seeds[0]).first);
            partner.entry.push_back(entries.at(seeds[1]).first);

            node   .box = entries.at(seeds[0]).second;
            partner.box = entries.at(seeds[1]).second;

            // remove them from entries pool
            entries.erase(entries.begin() + std::max(seeds[0], seeds[1]));
            entries.erase(entries.begin() + std::min(seeds[0], seeds[1]));
        }

        while(!entries.empty())
        {
            if(min_entry > node.entry.size() &&
               min_entry - node.entry.size() >= entries.size())
            {
                for(const auto& idx_box : entries)
                {
                    node.entry.push_back(idx_box.first);
                    node.box = this->expand(node.box, idx_box.second);
                }
                return partner;
            }
            if(min_entry > partner.entry.size() &&
               min_entry - partner.entry.size() >= entries.size())
            {
                for(const auto& idx_box : entries)
                {
                    partner.entry.push_back(idx_box.first);
                    partner.box = this->expand(partner.box, idx_box.second);
                }
                return partner;
            }

            const auto next = this->pick_next(entries, node.box, partner.box);
            if(next.second) // next is for node
            {
                node.entry.push_back(entries.at(next.first).first);
                node.box = this->expand(node.box, entries.at(next.first).second);
            }
            else // next is for partner
            {
                partner.entry.push_back(entries.at(next.first).first);
                partner.box = this->expand(partner.box, entries.at(next.first).second);
            }
            entries.erase(entries.begin() + next.first);
        }
        return partner;
    }

    // auxiliary function for the quadratic algorithm.
    std::array<std::size_t, 2>
    pick_seeds(const boost::container::static_vector<
            std::pair<std::size_t, box_type>, max_entry+1>& entries) const
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
    std::pair<std::size_t, bool> pick_next(const boost::container::static_vector<
            std::pair<std::size_t, box_type>, max_entry+1>& entries,
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
    boost::optional<std::pair<std::size_t,
        typename boost::container::static_vector<std::size_t, max_entry>::const_iterator>>
    find_leaf(std::size_t node_idx, const value_type& entry) const
    {
        const node_type& node = node_at(node_idx);
        const auto tight_box  = box_getter_(entry.second, 0.0);

        if(!(this->is_inside(tight_box, node.box)))
        {
            return boost::none;
        }

        if(node.is_leaf)
        {
            for(auto i=node.entry.begin(), e=node.entry.end(); i!=e; ++i)
            {
                // if the ID is the same, then the objects should be the same.
                if(container_.at(*i).first == entry.first)
                {
                    return std::make_pair(node_idx, i);
                }
            }
            return boost::none;
        }
        else // node is an internal node
        {
            for(auto i=node.entry.begin(), e=node.entry.end(); i!=e; ++i)
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
        assert(node.is_leaf);

        if(node.has_enough_entry() || node.parent == nil)
        {
            return; // nothing is needed.
        }

        // Leaf has too few entries. The tree structure should be balanced.
        // Here, the naivest approach is chosen. First remove the leaf node
        // that has too few children and re-insert all of them later.

        // copy objects in the leaf node that is being removed to re-insert them later
        std::vector<value_type> eliminated_objs;
        eliminated_objs.reserve(node.entry.size());
        for(const auto& idx: node.entry)
        {
            eliminated_objs.push_back(this->container_.at(idx));
            this->erase_value(idx);
        }

        const auto parent_idx = node.parent;

        // erase the node N from its parent and condense aabb
        auto found = std::find(this->node_at(node.parent).entry.begin(),
                               this->node_at(node.parent).entry.end(), N);
        assert(found != this->node_at(node.parent).entry.end());

        this->erase_node(*found);
        this->node_at(parent_idx).entry.erase(found);

        // condense node parent box without the leaf node N.
        this->condense_box(this->node_at(parent_idx));
        this->adjust_tree(parent_idx);

        // re-insert entries that were in node N
        for(const auto obj : eliminated_objs)
        {
            this->insert(obj);
        }
        // condense ancester nodes...
        this->condense_node(parent_idx);
        return;
    }

    // check the number of entries in the Nth internal node. If it has too few
    // entries, it removes the node and balance the tree.
    void condense_node(const std::size_t N)
    {
        const node_type& node = this->node_at(N);
        assert(node.is_leaf == false);

        if(node.has_enough_entry())
        {
            return; // if the node has enough number of entries, then it's okay.
        }
        if(node.parent == nil && node.entry.size() == 1)
        {
            // if the root has only one entry, then the child node of the current
            // root node should be the root node, no?
            this->root_ = node.entry.front();
            this->erase_node(N);
            assert(!this->is_valid_node_index(N));
            return;
        }

        // collect index of nodes that are children of the node to be removed
        const std::vector<std::size_t> eliminated_nodes(
                node.entry.begin(), node.entry.end());

        const auto parent_idx = node.parent;

        // erase the node N from its parent and condense its aabb
        auto found = std::find(this->node_at(parent_idx).entry.begin(),
                               this->node_at(parent_idx).entry.end(), N);
        assert(found != this->node_at(parent_idx).entry.end());

        // remove the node from its parent.entry
        this->erase_node(*found);
        this->node_at(parent_idx).entry.erase(found);
        this->condense_box(this->node_at(parent_idx));

        // re-insert nodes eliminated from node N
        for(const auto idx : eliminated_nodes)
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
        assert(!node.entry.empty());

        const auto& entries = node.entry;
        if(node.is_leaf)
        {
            node.box = box_getter_(this->container_.at(entries.front()).second,
                                   this->margin_);
            for(auto i = std::next(entries.begin()); i != entries.end(); ++i)
            {
                node.box = this->expand(node.box,
                    box_getter_(container_.at(*i).second, this->margin_));
            }
            return;
        }
        else
        {
            node.box = this->node_at(node.entry.front()).box;
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

        if(node_at(L).has_enough_storage())
        {
            node_at(L).entry.push_back(N);
            node_at(N).parent = L;
            node_at(L).box    = this->expand(node_at(L).box, entry);
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
            throw std::logic_error("the root has its parent node!?");
        }

        while(this->level_of(node_idx) != level)
        {
            Real diff_area_min = std::numeric_limits<Real>::max();
            Real area_min      = std::numeric_limits<Real>::max();

            const node_type& node = this->node_at(node_idx);
            for(const auto& entry_idx : node.entry)
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

    // the naivest approach; going down until it hits leaf
    std::size_t level_of(std::size_t node_idx) const
    {
        std::size_t level = 0;
        while(!(node_at(node_idx).is_leaf))
        {
            ++level;
            node_idx = node_at(node_idx).entry.front();
        }
        return level;
    }


private:

    // ========================================================================
    // utility member methods to handle `container_` and `tree_`.
    //
    // RTree represents nodes and values by their index in a container_.
    // Thus the indices should be kept when we insert a new node/value or
    // remove an old node/value.

    // ------------------------------------------------------------------------
    // for nodes

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
        node_at(i) = node_type(false, nil);
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

    // the same thing as above + assign it to rmap.
    std::size_t add_value(const value_type& v)
    {
        if(overwritable_values_.empty())
        {
            const std::size_t new_index = container_.size();
            container_.push_back(v);
            rmap_[v.first] = new_index;
            return new_index;
        }
        const std::size_t new_index = overwritable_values_.back();
        overwritable_values_.pop_back();
        container_.at(new_index) = v;
        rmap_[v.first] = new_index;
        return new_index;
    }
    // the same thing as above + remove corresponding key from rmap.
    void erase_value(const std::size_t i)
    {
        using std::swap;
        overwritable_values_.push_back(i);
        value_type old; // default value
        swap(this->container_.at(i), old);

        const auto status = rmap_.erase(old.first);
        assert(status == 1);
        (void)status;
        return ;
    }
    bool is_valid_value_index(const std::size_t i) const noexcept
    {
        return std::find(overwritable_values_.begin(),
                         overwritable_values_.end(), i) ==
            overwritable_values_.end();
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

        return box_type(lw + d, up + d);
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

    bool is_inside_of_boundary(const Real3 r) const noexcept
    {
        if(r[0] < 0 || edge_lengths_[0] < r[0]) {return false;}
        if(r[1] < 0 || edge_lengths_[1] < r[1]) {return false;}
        if(r[2] < 0 || edge_lengths_[2] < r[2]) {return false;}
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

    std::size_t           root_;                // index of the root node in tree_
    Real                  margin_;              // margin of the Box.
    Real3                 edge_lengths_;        // (Lx, Ly, Lz) of PBC
    tree_type             tree_;                // vector of nodes.
    container_type        container_;           // vector of values to contain.
    key_to_value_map_type rmap_;                // map: ObjectID -> Object
    index_buffer_type     overwritable_values_; // list of "already-removed" values
    index_buffer_type     overwritable_nodes_;  // ditto.
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
