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
#include <set>

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
    typedef AABB                                      box_type;
    typedef AABBGetter                                box_getter_type;
    typedef std::pair<box_type, std::size_t>          rtree_value_type;
    typedef std::pair<ObjectID, Object>               value_type;
    typedef std::vector<value_type>                   container_type;
    typedef std::unordered_map<ObjectID, std::size_t> key_to_value_map_type;

    static constexpr std::size_t min_entry = MinEntry;
    static constexpr std::size_t max_entry = MaxEntry;
    static constexpr std::size_t nil = std::numeric_limits<std::size_t>::max();

    struct node_type
    {
        typedef boost::container::static_vector<std::size_t, max_entry>
                entry_type;

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
    //     bool operator(const value_type&, const Real3& edge_lengths);
    //
    //     // to skip non-interacting nodes purely geometric criteria.
    //     bool operator(const AABB&,       const Real3& edge_lengths);
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

    bool update(const ObjectID& id, const Object& obj)
    {
        const auto v = std::make_pair(id, obj);
        if(!this->has(id))
        {
            this->insert(v);
            return true;
        }
        const auto old_v = this->container_.at(this->rmap_.at(id));

        const auto found = this->find_leaf(this->root_, old_v);
        if(!found)
        {
            throw std::logic_error("[error] internal error in PeriodicRTree");
        }
        const auto node_idx  =   found->first;
        const auto value_idx = *(found->second);

        const auto tight_box = this->box_getter_(obj, 0.0);
        const auto& node_box = this->tree_.at(node_idx).box;

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

            // erasing the object in a node.
            // remove the entry, then condence the AABB.
            this->tree_.at(node_idx).entry.erase(found->second);
            this->erase_value(value_idx);

            this->condense_box(this->tree_.at(node_idx));
            this->condense_leaf(node_idx);
            // done.

            this->insert(v);
        }
        return false;
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

        // TODO bug here maybe
        if(tree_.at(L).has_enough_storage())
        {
            tree_.at(L).entry.push_back(idx);
            tree_.at(L).box = this->expand(tree_.at(L).box, box);
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

    void erase(const value_type& v)
    {
        if(this->root_ == nil)
        {
            throw std::out_of_range("PeriodicRTree::erase: tree is empty.");
        }
        // find_leaf returns pair{leaf_node_idx, iterator of leaf_node.entry}
        if(const auto found = this->find_leaf(this->root_, v))
        {
            const auto node_idx  =   found->first;
            const auto value_idx = *(found->second);
            this->tree_.at(node_idx).entry.erase(found->second);
            this->erase_value(value_idx);

            this->condense_box(this->tree_.at(node_idx));
            this->condense_leaf(node_idx);
            return;
        }
        else
        {
            throw std::out_of_range("PeriodicRTree::erase: value not found.");
        }
    }

    Real3 const& edge_lengths() const noexcept {return edge_lengths_;}

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

    // user-defined AABBGetter may be stateful (unlikely).
    AABBGetter const& get_aabb_getter() const noexcept {return box_getter_;}

    // box margin.
    Real margin() const noexcept {return margin_;}

    bool diagnosis() const
    {
#define ECELL4_PERIODIC_RTREE_ASSERT(x)\
        do{\
            if(!(x)){this->dump(std::cerr);}\
            assert((x));\
        }while(false);

        std::size_t num_objects = 0;
        std::size_t num_inodes  = 1;
        for(std::size_t i=0; i<tree_.size(); ++i)
        {
            if(!this->is_valid_node_index(i)){continue;}

            if(this->tree_.at(i).is_leaf)
            {
                num_objects += this->tree_.at(i).entry.size();
            }
            else
            {
                num_inodes += this->tree_.at(i).entry.size();
            }
        }
        ECELL4_PERIODIC_RTREE_ASSERT(this->list_objects().size() == num_objects);
        ECELL4_PERIODIC_RTREE_ASSERT(this->tree_.size() - this->overwritable_nodes_.size() == num_inodes);

        bool root_found = false;
        for(std::size_t i=0; i<tree_.size(); ++i)
        {
            if(!this->is_valid_node_index(i)){continue;}

            if(this->tree_.at(i).parent == nil)
            {
                ECELL4_PERIODIC_RTREE_ASSERT(!root_found);
                root_found = true;
            }
            else
            {
                const auto& e = this->tree_.at(this->tree_.at(i).parent).entry;
                ECELL4_PERIODIC_RTREE_ASSERT(this->is_valid_node_index(this->tree_.at(i).parent));
                ECELL4_PERIODIC_RTREE_ASSERT(!this->tree_.at(this->tree_.at(i).parent).is_leaf);
                if(std::find(e.begin(), e.end(), i) == e.end())
                {
                    std::cerr << "node " << i << " is not found in the list of entries in parent, " << this->tree_.at(i).parent << std::endl;
                }
                ECELL4_PERIODIC_RTREE_ASSERT(std::find(e.begin(), e.end(), i) != e.end());
            }
        }

        for(std::size_t i=0; i<tree_.size(); ++i)
        {
            if(!this->is_valid_node_index(i)) {continue;}
            if(tree_.at(i).is_leaf && !diagnosis_rec(i))
            {
                std::cerr << "AABB does not wrap the leaf node" << std::endl;
                return false;
            }
        }

        for(std::size_t i=0; i<container_.size(); ++i)
        {
            if(!this->is_valid_value_index(i)) {continue;}

            if(rmap_.count(container_.at(i).first) != 0)
            {
                if(rmap_.at(container_.at(i).first) != i)
                {
                    std::cerr << "key->index map broken: should be " << i << " but got " << rmap_.at(container_.at(i).first) << std::endl;
                    return false;
                }
            }

            if(!this->find_leaf(this->root_, container_.at(i)))
            {
                std::cerr << "object with ID " << container_.at(i).first << " is not found by find_leaf" << std::endl;
                return false;
            }
        }
#undef ECELL4_PERIODIC_RTREE_ASSERT
        return true;
    }

    bool diagnosis_rec(std::size_t node_idx) const
    {
        while(tree_.at(node_idx).parent != nil)
        {
            const auto& node = this->tree_.at(node_idx);
            if(!(this->is_inside(node.box, tree_.at(node.parent).box)))
            {
                std::cerr << "AABB of node " << node_idx << ", " << node.box
                          << ", is not within that of its parent, " << tree_.at(node.parent).box << "." << std::endl;
                return false;
            }
            node_idx = node.parent;
        }
        return true;
    }

    void dump(std::ostream& os) const
    {
        os << "-------------------------------------------------------\n";
        std::ostringstream oss;
        if(this->is_valid_node_index(this->root_)) {oss << ' ';} else {oss << '!';}
        oss << '(' << std::setw(3) << this->root_ << ") ";
        this->dump_rec(this->root_, os, oss.str());
        os << "-------------------------------------------------------\n";
        os << std::flush;
        return ;
    }

private:

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
    // recursively check nodes. if the node is a leaf, check values inside it.

    template<typename F, typename OutputIterator>
    OutputIterator query_recursive(
            const std::size_t node_idx, F matches, OutputIterator& out) const
    {
        const node_type& node = this->tree_.at(node_idx);
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
            const auto& node_aabb = tree_.at(entry).box;
            if(matches(node_aabb, this->edge_lengths_))
            {
                this->query_recursive(entry, matches, out);
            }
        }
        return out;
    }

private:

    // ------------------------------------------------------------------------
    // construct and manage RTree structure.

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
        while(!(this->tree_.at(node_idx).is_leaf))
        {
            // find minimum expansion
            Real diff_area_min = std::numeric_limits<Real>::max();
            Real area_min      = std::numeric_limits<Real>::max();

            const auto& node = this->tree_.at(node_idx);
            for(const std::size_t i : node.entry)
            {
                const auto& current_box = this->tree_.at(i).box;

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

    // node has another value.
    void adjust_tree(std::size_t node_idx)
    {
        while(tree_.at(node_idx).parent != nil)
        {
            const auto& node   = this->tree_.at(node_idx);
            auto&       parent = this->tree_.at(node.parent);

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
    void adjust_tree(const std::size_t N, const std::size_t NN)
    {
        assert(N != NN);
        // we hit the root. to assign a new node, we need to make tree deeper.
        if(tree_.at(N).parent == nil)
        {
            node_type new_root(/*is_leaf = */false, /*parent = */ nil);
            new_root.entry.push_back(N);
            new_root.entry.push_back(NN);
            new_root.box = this->expand(tree_.at(N).box, tree_.at(NN).box);
            this->root_  = this->add_node(std::move(new_root));

            this->tree_.at( N).parent = this->root_;
            this->tree_.at(NN).parent = this->root_;
            return;
        }
        else
        {
            const auto& node    = tree_.at(N);
            const auto& partner = tree_.at(NN);
            assert(node.parent == partner.parent);

            const auto parent_idx = node.parent;
            auto& parent = tree_.at(parent_idx);
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

    // split nodes by quadratic algorithm
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

        const std::size_t PP = this->add_node(node_type(false, tree_.at(P).parent));
        assert(P  != PP);
        assert(NN != PP);
        node_type& node    = tree_.at(P);
        node_type& partner = tree_.at(PP);

        assert(!node.is_leaf);
        assert(!partner.is_leaf);

        boost::container::static_vector<
            std::pair<std::size_t, box_type>, max_entry + 1> entries;
        entries.emplace_back(NN, tree_.at(NN).box);

        for(const auto& entry_idx : node.entry)
        {
            entries.emplace_back(entry_idx, tree_.at(entry_idx).box);
        }
        node.entry.clear();
        partner.entry.clear();

        // assign first 2 entries to node and partner
        {
            const auto seeds = this->pick_seeds(entries);
            assert(seeds[0] != seeds[1]);

            node   .entry.push_back(entries.at(seeds[0]).first);
            partner.entry.push_back(entries.at(seeds[1]).first);

            tree_.at(entries.at(seeds[0]).first).parent = P;
            tree_.at(entries.at(seeds[1]).first).parent = PP;

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
                    tree_.at(idx_box.first).parent = P;
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
                    tree_.at(idx_box.first).parent = PP;
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
                tree_.at(entries.at(next.first).first).parent = P;
                node.box = this->expand(node.box, entries.at(next.first).second);
            }
            else
            {
                partner.entry.push_back(entries.at(next.first).first);
                tree_.at(entries.at(next.first).first).parent = PP;
                partner.box = this->expand(partner.box, entries.at(next.first).second);
            }
            entries.erase(entries.begin() + next.first);
        }
        tree_.at(P)  = node;
        tree_.at(PP) = partner;
        return PP;
    }

    node_type split_leaf(const std::size_t N, const std::size_t vidx,
                         const box_type& entry)
    {
        node_type& node = tree_.at(N);
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
        const node_type& node = tree_.at(node_idx);
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
                if(!(this->is_inside(tight_box, this->tree_.at(*i).box)))
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

    void condense_leaf(const std::size_t N)
    {
        const node_type& node = this->tree_.at(N);
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
        assert(this->is_valid_node_index(parent_idx));

        // erase the node N from its parent and condense aabb
        auto found = std::find(this->tree_.at(node.parent).entry.begin(),
                               this->tree_.at(node.parent).entry.end(), N);
        assert(found != this->tree_.at(node.parent).entry.end());

        this->erase_node(*found);
        assert(!this->is_valid_node_index(*found));
        this->tree_.at(parent_idx).entry.erase(found);

        // condense node parent box without the leaf node N.
        this->condense_box(this->tree_.at(parent_idx));

        // re-insert entries that were in node N
        for(const auto obj : eliminated_objs)
        {
            this->insert(obj);
        }

        // condense ancester nodes...
        this->condense_node(parent_idx);
        return;
    }

    void condense_node(const std::size_t N)
    {
        const node_type& node = this->tree_.at(N);
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
        auto found = std::find(this->tree_.at(parent_idx).entry.begin(),
                               this->tree_.at(parent_idx).entry.end(), N);
        assert(found != this->tree_.at(parent_idx).entry.end());

        // remove the node from its parent.entry
        this->erase_node(*found);
        this->tree_.at(parent_idx).entry.erase(found);
        this->condense_box(this->tree_.at(parent_idx));

        // re-insert nodes eliminated from node N
        for(const auto idx : eliminated_nodes)
        {
            this->re_insert(idx);
        }
        this->condense_node(parent_idx);
        return;
    }

    // re-calculate the AABB of a node
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
            node.box = this->tree_.at(node.entry.front()).box;
            for(auto i = std::next(entries.begin()); i != entries.end(); ++i)
            {
                node.box = this->expand(node.box, this->tree_.at(*i).box);
            }
        }
        return;
    }

    void re_insert(const std::size_t N)
    {
        // reset connection to the parent!

        // insert node to its proper parent. to find the parent of this node N,
        // add 1 to level. root node should NOT come here.
        const std::size_t level = this->level_of(N) + 1;
        const box_type&   entry = this->tree_.at(N).box;
        const std::size_t L = this->choose_node_with_level(entry, level);

        if(tree_.at(L).has_enough_storage())
        {
            tree_.at(L).entry.push_back(N);
            tree_.at(N).parent = L;
            tree_.at(L).box    = this->expand(tree_.at(L).box, entry);
            this->adjust_tree(L);
        }
        else
        {
            const std::size_t LL = this->split_node(L, N);
            this->adjust_tree(L, LL);
        }
        return;
    }

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

            const node_type& node = this->tree_.at(node_idx);
            for(const auto& entry_idx : node.entry)
            {
                const auto& entry_box = tree_.at(entry_idx).box;
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

    // the naive approach; going down until it hits leaf
    std::size_t level_of(std::size_t node_idx) const
    {
        std::size_t level = 0;
        while(!(tree_.at(node_idx).is_leaf))
        {
            ++level;
            node_idx = tree_.at(node_idx).entry.front();
        }
        return level;
    }


private:

    // ------------------------------------------------------------------------
    // utility member methods to handle `container_` and `tree_`.
    //
    // RTree represents nodes and values by their index in a container_.
    // Thus the indices should be kept when we insert a new node/value or
    // remove an old node/value.

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
        tree_.at(new_index) = n;
        return new_index;
    }
    // mark index `i` overritable and fill old value by the default value.
    void erase_node(const std::size_t i)
    {
        overwritable_nodes_.push_back(i);
        tree_.at(i) = node_type(false, nil);
        return;
    }

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

    bool is_valid_node_index(const std::size_t i) const noexcept
    {
        return std::find(overwritable_nodes_.begin(),
                         overwritable_nodes_.end(), i) ==
            overwritable_nodes_.end();
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

        assert(this->is_inside_of_boundary(lc));
        assert(this->is_inside_of_boundary(rc));
        assert(lhs.lower()[0] <= lhs.upper()[0]);
        assert(lhs.lower()[1] <= lhs.upper()[1]);
        assert(lhs.lower()[2] <= lhs.upper()[2]);
        assert(rhs.lower()[0] <= rhs.upper()[0]);
        assert(rhs.lower()[1] <= rhs.upper()[1]);
        assert(rhs.lower()[2] <= rhs.upper()[2]);

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
        const auto dc = abs(this->restrict_direction(lc - rc));

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
