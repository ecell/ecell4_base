#ifndef VMAP_HPP
#define VMAP_HPP 

#include <cstring>
#include <utility>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/type_traits.hpp>
#include <boost/call_traits.hpp>
#include <boost/mpl/if.hpp>

#include "utils.hpp"

template<typename Tkey_, typename Tval_,
        template<typename, typename> class MFget_mapper_ =
            get_default_impl::std::template map,
        template<typename> class MFget_racntnr_ =
            get_default_impl::std::template vector>
class vmap
{
public:
    template<typename Tfirst_, typename Tsecond_>
    struct pair
    {
        typedef Tfirst_ first_type;
        typedef Tsecond_ second_type;

        pair(typename boost::call_traits<Tfirst_>::param_type _first,
                typename boost::call_traits<Tsecond_>::param_type _second)
                : first(_first), second(_second) {}

        Tfirst_ first;
        Tsecond_ second;
    };

    typedef typename MFget_racntnr_<Tval_>::type
            random_accessible_container_type;
    typedef typename MFget_mapper_<Tkey_,
            typename random_accessible_container_type::size_type>::type
                 mapper_type;
    typedef typename MFget_racntnr_<typename mapper_type::key_type>::type
            reverse_mapper_type;
    typedef Tval_ mapped_type;
    typedef Tkey_ key_type;
    typedef typename mapper_type::size_type size_type;
    typedef typename mapper_type::difference_type difference_type;
    typedef std::pair<const Tkey_, Tval_> value_type;
    typedef pair<const Tkey_&, Tval_&> reference;
    typedef pair<const Tkey_&, const Tval_&> const_reference;
    typedef typename random_accessible_container_type::iterator value_iterator;
    typedef typename random_accessible_container_type::const_iterator const_value_iterator;
    typedef boost::iterator_range<value_iterator> value_range;
    typedef boost::iterator_range<const_value_iterator> const_value_range;

    template<typename Thost_, typename Treftype_, typename Tcreftype_,
                typename Tmiter_, typename Tracntnr_>
    class iterator_base
    {
        friend class vmap;
    public:
        typedef Treftype_ reference;
        typedef Tcreftype_ const_reference;

    public:
        iterator_base(const Tmiter_& miter,
                Tracntnr_& racntnr)
                : miter_(miter), racntnr_(racntnr) {}

        iterator_base(const iterator_base& that)
                : miter_(that.miter_), racntnr_(that.racntnr_) {}

        bool operator==(const Thost_& rhs)
        {
            return miter_ == rhs.miter_;
        }

        bool operator!=(const Thost_& rhs)
        {
            return miter_ != rhs.miter_;
        }

        friend Thost_& operator++(Thost_& self)
        {
            ++self.miter_;
            return self;
        }

        friend Thost_& operator--(Thost_& self)
        {
            --self.miter_;
            return self;
        }

        friend Thost_ operator++(Thost_& self, int)
        {
            Thost_ retval(self);
            ++self.miter_;
            return retval;
        }

        friend Thost_ operator--(Thost_& self, int)
        {
            Thost_ retval(self);
            --self.miter_;
            return retval;
        }

        reference operator*()
        {
            return reference((*miter_).first, racntnr_[(*miter_).second]);
        }

        const_reference operator*() const
        {
            return const_reference((*miter_).first, racntnr_[(*miter_).second]);
        }

        Thost_& operator=(const iterator_base& rhs)
        {
            // XXX: this works ;-p
            std::memmove(this, &rhs, sizeof(rhs));
            return *reinterpret_cast<Thost_*>(this);
        }

    public:
        Tmiter_ miter_;
        Tracntnr_& racntnr_;
    };

    class iterator: public iterator_base<iterator, reference, const_reference,
            typename mapper_type::iterator,
            random_accessible_container_type>
    {
    private:
        typedef iterator_base<iterator, reference, const_reference,
                typename mapper_type::iterator,
                random_accessible_container_type> base_type;
    public:
        iterator(const typename mapper_type::iterator& miter,
                random_accessible_container_type& racntnr)
            : base_type(miter, racntnr) {}

        iterator(const iterator& that): base_type(that) {}
    };

    class const_iterator
        : public iterator_base<const_iterator, const_reference, const_reference,
            typename mapper_type::const_iterator,
            const random_accessible_container_type>
    {
    private:
        typedef iterator_base<const_iterator, const_reference, const_reference,
                typename mapper_type::const_iterator,
                const random_accessible_container_type> base_type;

    public:
        const_iterator(const typename mapper_type::const_iterator& miter,
                const random_accessible_container_type& racntnr)
            : base_type(miter, racntnr) {}

        const_iterator(const const_iterator& that)
            : base_type(that) {}

        const_iterator(const iterator& that)
            : base_type(that.miter_, that.racntnr_) {}
    };

public:
    vmap() {}

    size_type erase(const key_type& k)
    {
        return erase(mapper_.equal_range(k));
    }

    void erase(const iterator& p, const iterator& q)
    {
        erase(std::make_pair(p.miter_, q.miter_));
    }

    void erase(const iterator& p)
    {
        erase(std::make_pair(p.miter_,
                ++typename mapper_type::iterator(p.miter_)));
    }

    void clear()
    {
        racntnr_.clear();
        rmapper_.clear();
        mapper_.clear();
    }

    std::pair<const_iterator, const_iterator> equal_range(const key_type& k) const
    {
        std::pair<typename mapper_type::const_iterator,
                typename mapper_type::const_iterator> r(mapper_.equal_range(k));
        return std::make_pair(
                const_iterator(r.begin(), racntnr_),
                const_iterator(r.end(), racntnr_));
    }

    std::pair<iterator, iterator> equal_range(const key_type& k)
    {
        std::pair<typename mapper_type::const_iterator,
                typename mapper_type::const_iterator> r(mapper_.equal_range(k));
        return std::make_pair(
                iterator(r.begin(), racntnr_),
                iterator(r.end(), racntnr_));
    }

    const_iterator find(const key_type& k) const
    {
        return const_iterator(mapper_.find(k), racntnr_);
    }

    iterator find(const key_type& k)
    {
        return iterator(mapper_.find(k), racntnr_);
    }

    const_iterator begin() const
    {
        return const_iterator(mapper_.begin(), racntnr_);
    }

    iterator begin()
    {
        return iterator(mapper_.begin(), racntnr_);
    }

    const_iterator end() const
    {
        return const_iterator(mapper_.end(), racntnr_);
    }

    iterator end()
    {
        return iterator(mapper_.end(), racntnr_);
    }

    size_type count(const key_type& k)
    {
        return mapper_.count(k);
    }

    value_range values()
    {
        return value_range(racntnr_.begin(), racntnr_.end());
    }

    const_value_range values() const
    {
        return value_range(racntnr_.begin(), racntnr_.end());
    }

    std::pair<iterator, bool> insert(const value_type& v)
    {
        std::pair<typename mapper_type::iterator, bool> ir(
                mapper_.insert(
                    std::make_pair(v.first, racntnr_.size())));
        if (!ir.second)
        {
            racntnr_[(*ir.first).second] = v.second;
            return make_pair(iterator(ir.first, racntnr_), false);
        }

        racntnr_.push_back(v.second);
        rmapper_.push_back(v.first);
        return make_pair(iterator(ir.first, racntnr_), true);
    }

    size_type size() const
    {
        return mapper_.size();
    }

    mapped_type& operator[](const key_type& key)
    {
        std::pair<typename mapper_type::iterator, bool> ir(
                mapper_.insert(std::make_pair(key, racntnr_.size())));
        if (!ir.second) {
            return racntnr_[(*ir.first).second];
        }

        racntnr_.push_back(mapped_type());
        rmapper_.push_back(key);
        return racntnr_.back();
    }

private:
    size_type erase(const std::pair<typename mapper_type::iterator,
                typename mapper_type::iterator>& r)
    {
        typename random_accessible_container_type::size_type b(racntnr_.size());

        for (typename mapper_type::iterator i(r.first);
                i != r.second; ++i)
        {
            --b;
            const key_type& k(rmapper_[b]);
            std::swap(racntnr_[b], racntnr_[(*i).second]);
            std::swap(rmapper_[b], rmapper_[(*i).second]);
            std::swap(mapper_[k], (*i).second);
            mapper_.erase((*i).first); 
        }
        racntnr_.resize(b);
        rmapper_.resize(b);
        return racntnr_.size() - b;
    }

private:
    mapper_type mapper_;
    reverse_mapper_type rmapper_;
    random_accessible_container_type racntnr_;
};

#endif /* VMAP_HPP */
