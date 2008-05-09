#ifndef _OBJECT_CONTAINER_HPP
#define _OBJECT_CONTAINER_HPP

#include <tr1/unordered_map>
#include <boost/coroutine/generator.hpp>
#include <boost/bind.hpp>
#include "object_container.hpp"
#include "filters.hpp"

typedef object_container<double, int, make_get_mapper_mf<std::tr1::unordered_map>::meta_type> my_object_container_type;

struct sphere_ref
{
    typedef my_object_container_type::iterator impl_type;

public:
    sphere_ref(impl_type impl): impl_(impl) {}

    const double& x() const
    {
        return (*impl_).second.position.x();
    }

    double& x()
    {
        return (*impl_).second.position.x();
    }

    const double& y() const
    {
        return (*impl_).second.position.y();
    }

    double& y()
    {
        return (*impl_).second.position.y();
    }

    const double& z() const
    {
        return (*impl_).second.position.z();
    }

    double& z()
    {
        return (*impl_).second.position.z();
    }

    const double& radius() const
    {
        return (*impl_).second.radius;
    }

    double& radius()
    {
        return (*impl_).second.radius;
    }

    int id() const
    {
        return (*impl_).first;
    }

    operator sphere<double>() const
    {
        return (*impl_).second;
    }

private:
    impl_type impl_;
};

template<typename Tstrm_>
inline std::basic_ostream<Tstrm_>& operator<<(std::basic_ostream<Tstrm_>& strm,
        const sphere_ref& v)
{
    strm << static_cast<sphere<double> >(v);
    return strm;
}

sphere<double> *__impl_sphere_ref_clone(const sphere_ref* src)
{
    return new sphere<double>(*src);
}

struct take_neighbor_collector
{
    typedef std::pair<sphere_ref*,
            my_object_container_type::position_type::value_type> value_type;

    typedef boost::coroutines::generator<const value_type*> generator_type;

    inline take_neighbor_collector(generator_type::self& self)
            : self_(self)
    {
    }

    inline void operator()(my_object_container_type::iterator i,
            const my_object_container_type::position_type::value_type& d)
    {
        last_.first = new sphere_ref(i);
        last_.second = d;
        self_.yield(&last_);
    }

private:
    generator_type::self& self_;
    value_type last_;
};

inline const take_neighbor_collector::value_type*
take_neighbor_generator(
        take_neighbor_collector::generator_type::self& self,
        my_object_container_type& cntnr,
        const my_object_container_type::mapped_type& sphere)
{
    take_neighbor_collector col(self);
    take_neighbor(cntnr, col, sphere);
    self.yield(NULL);
    self.exit();
    return NULL; // never get here
}

inline const take_neighbor_collector::value_type*
take_neighbor_cyclic_generator(
        take_neighbor_collector::generator_type::self& self,
        my_object_container_type& cntnr,
        const my_object_container_type::mapped_type& sphere)
{
    take_neighbor_collector col(self);
    take_neighbor_cyclic(cntnr, col, sphere);
    self.yield(NULL);
    self.exit();
    return NULL; // never get here
}

inline take_neighbor_collector::generator_type*
__impl_object_container_iterneighbors(my_object_container_type* pimpl,
        const my_object_container_type::mapped_type* psphere)
{
    return new take_neighbor_collector::generator_type(
            boost::bind(take_neighbor_generator, _1, *pimpl, *psphere));
}

inline take_neighbor_collector::generator_type*
__impl_object_container_iterneighbors_cyclic(my_object_container_type* pimpl,
        const my_object_container_type::mapped_type* psphere)
{
    return new take_neighbor_collector::generator_type(
            boost::bind(take_neighbor_cyclic_generator, _1, *pimpl, *psphere));
}

inline sphere_ref* __impl_object_container_find(
        my_object_container_type* pimpl,
        const my_object_container_type::key_type& key)
{
    my_object_container_type::iterator i(pimpl->find(key));
    if (i == pimpl->end())
    {
        return NULL;
    }
    return new sphere_ref(i);
}

inline bool
__impl_object_container_insert(my_object_container_type* cntnr,
        const my_object_container_type::key_type& key,
        my_object_container_type::mapped_type* sphere)
{
    return cntnr->insert(my_object_container_type::value_type(key, *sphere)).second;
}


#endif /* _OBJECT_CONTAINER_HPP */
