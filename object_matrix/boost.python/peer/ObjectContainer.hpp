#ifndef OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP
#define OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP

#include <functional>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/add_const.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <numpy/arrayobject.h>

#include "peer/utils.hpp"
#include "peer/numpy/type_mappings.hpp"

#include "peer/tuple_converters.hpp"
#include "peer/numpy/pyarray_backed_allocator.hpp"
#include "peer/numpy/ndarray_converters.hpp"

#include "position.hpp"
#include "sphere.hpp"
#include "filters.hpp"

#include "get_mapper_mf.hpp"

#if OBJECTMATRIX_USE_ITERATOR
#include <boost/coroutine/generator.hpp>
#include <boost/bind.hpp>
#include "peer/generator_support.hpp"
#endif /* OBJECTMATRIX_USE_ITERATOR */

template<typename Tval_>
struct get_mapper_mf<boost::python::object, Tval_>
{
#if HAVE_UNORDERED_MAP || HAVE_TR1_UNORDERED_MAP || HAVE_EXT_HASH_MAP
    struct hasher: public std::unary_function<boost::python::object, std::size_t>
    {
        typedef boost::python::object argument_type;
        typedef std::size_t result_type;

        result_type operator()(const argument_type& arg) const
        {
            return static_cast<result_type>((long)PyObject_Hash(arg.ptr()));
        }
    };
#endif

#if HAVE_UNORDERED_MAP
    typedef std::unordered_map<boost::python::object, Tval_, hasher> type;
#elif HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<boost::python::object, Tval_, hasher> type;
#elif HAVE_EXT_HASH_MAP
    typedef __gnu_cxx::hash_map<boost::python::object, Tval_, hasher> type;
#else 
    typedef std::map<boost::python::object, Tval_> type;
#endif
};

namespace peer {

struct ObjectContainerBase
{
    template<typename T_, typename Tref_>
    struct add_const_if:
        boost::mpl::if_<boost::is_const<Tref_>,
            typename boost::add_const<T_>::type, T_> {};

    template<typename Tpair_>
    struct select1st
        : public std::unary_function<
                Tpair_&,
                typename add_const_if<
                    typename Tpair_::first_type, Tpair_>::type&>
    {
        typename add_const_if<
            typename Tpair_::first_type, Tpair_>::type&
        operator()(Tpair_& pair) const
        {
            return pair.first;
        }
    };
};

template<typename Timpl_>
class ObjectContainer: public ObjectContainerBase
{
public:
    typedef Timpl_ impl_type;
    typedef typename impl_type::key_type key_type;
    typedef typename impl_type::mapped_type mapped_type;
    typedef typename impl_type::position_type position_type;
    typedef typename impl_type::length_type length_type;
    typedef typename impl_type::size_type size_type;
    typedef typename impl_type::matrix_type::size_type matrix_size_type;

#ifdef OBJECTMATRIX_USE_ITERATOR
    class Generators
    {
    public:
        typedef std::pair<impl_type::iterator, position_type::value_type>
                result_type;

        typedef boost::coroutines::generator<const result_type*> generator_type; 
    private:
        struct collector: public std::binary_function<
                impl_type::reference, position_type::value_type, void>
        {
        public:
            typedef impl_type::iterator first_argument_type;
            typedef position_type::value_type second_argument_type;
            typedef void result_type;

        public:
            inline collector(generator_type::self& self,
                    impl_type::iterator end)
                    : self_(self), last_(end, 0)
            {
            }

            inline void operator()(impl_type::iterator i,
                    const position_type::value_type& d)
            {
                last_.first = i;
                last_.second = d;
                self_.yield(&last_);
            }
private:
            generator_type::self& self_;
            Generators::result_type last_;
        };

    public:
        inline static Enumerator<const result_type&>* enumerate_neighbors(
            ObjectContainer& cntnr, const mapped_type& sphere)
        {
            return util::make_enumerator(generator_type(
                boost::bind(&gen_take_neighbor, _1, cntnr, sphere)));
        }

        inline static Enumerator<const result_type&>* enumerate_neighbors_cyclic(
            ObjectContainer& cntnr, const mapped_type& sphere)
        {
            return util::make_enumerator(generator_type(
                boost::bind(&gen_take_neighbor_cyclic, _1, cntnr, sphere)));
        }

    private:
        Generators() {}

        inline static const result_type* gen_take_neighbor(
                generator_type::self& self,
                ObjectContainer& cntnr,
                const mapped_type& sphere)
        {
            collector col(self, cntnr.impl_.end());
            ::take_neighbor(cntnr.impl_, col, sphere);
            self.yield(NULL);
            self.exit();
            return NULL; // never get here
        }

        inline static const result_type* gen_take_neighbor_cyclic(
                generator_type::self& self,
                ObjectContainer& cntnr,
                const mapped_type& sphere)
        {
            collector col(self, cntnr.impl_.end());
            ::take_neighbor_cyclic(cntnr.impl_, col, sphere);
            self.yield(NULL);
            self.exit();
            return NULL; // never get here
        }
    };
#endif /* OBJECTMATRIX_USE_ITERATOR */

    class Builders
    {
    public:
        typedef std::vector<typename impl_type::key_type> key_array_type;
        typedef std::vector<length_type, util::pyarray_backed_allocator<length_type> >
                distance_array_type;
        typedef boost::tuple<key_array_type, distance_array_type>
                result_type;

        struct collector: public std::binary_function<
                typename impl_type::reference,
                typename position_type::value_type, void>
        {
            typedef typename impl_type::iterator first_argument_type;
            typedef typename position_type::value_type second_argument_type;
            typedef void result_type;
        public:
            inline collector(typename Builders::result_type& result)
                : //sa_(boost::get<0>(result)),
                ka_(boost::get<0>(result)),
                da_(boost::get<1>(result)) {}

            inline void operator()(typename impl_type::iterator i,
                    const typename position_type::value_type& d)
            {
                //sa_.push_back(i);
                ka_.push_back((*i).first);
                da_.push_back(d);
            }

        private:
            key_array_type& ka_;
            distance_array_type& da_;
        };

        struct all_neighbors_collector: public std::binary_function<
                typename impl_type::reference, typename position_type::value_type, void>
        {
            typedef typename impl_type::iterator first_argument_type;
            typedef typename position_type::value_type second_argument_type;
            typedef void result_type;
        public:
            inline all_neighbors_collector(typename Builders::result_type& result,
                    const position_type& pos)
                : //sa_(boost::get<0>(result)),
                ka_(boost::get<0>(result)),
                      da_(boost::get<1>(result)),
                      pos_(pos) {}

            inline void operator()(typename impl_type::iterator i)
            {
                //sa_.push_back(i);
                ka_.push_back((*i).first);
                da_.push_back(pos_.distance((*i).second.position()) 
                              - (*i).second.radius());

            }

            inline void operator()(typename impl_type::iterator i,
                    const position_type& d)
            {
                //sa_.push_back(i);
                ka_.push_back((*i).first);
                da_.push_back(pos_.distance((*i).second.position() + d)
                              - (*i).second.radius());
            }

        private:
            key_array_type& ka_;
            distance_array_type& da_;
            position_type pos_;
        };


    public:
        inline static void
        build_neighbors_array(result_type& retval,
                              ObjectContainer& cntnr,
                              const ::sphere<length_type>& sphere)
        {
            collector col(retval);
            take_neighbor(cntnr.impl_, col, sphere);
        }

        inline static void
        build_neighbors_array_cyclic(result_type& retval,
                ObjectContainer& cntnr, const ::sphere<length_type>& sphere)
        {
            collector col(retval);
            take_neighbor_cyclic(cntnr.impl_, col, sphere);
        }

        inline static void
        build_all_neighbors_array(result_type& retval,
                ObjectContainer& cntnr, const position_type& pos)
        {
            all_neighbors_collector col(retval, pos);
            cntnr.impl_.each_neighbor(cntnr.impl_.index(pos), col);
        }

        inline static void
        build_all_neighbors_array_cyclic(result_type& retval,
                ObjectContainer& cntnr, const position_type& pos)
        {
            all_neighbors_collector col(retval, pos);
            cntnr.impl_.each_neighbor_cyclic(cntnr.impl_.index(pos), col);
        }

    private:
        Builders() {}
    };

public:
    ObjectContainer(typename impl_type::length_type world_size,
            typename impl_type::matrix_type::size_type size)
            : impl_(world_size, size) {}

    size_type __len__() const
    {
        return impl_.size();
    }

    size_type getMatrixSize() const
    {
        return impl_.matrix_size();
    }

    length_type getWorldSize() const
    {
        return impl_.world_size();
    }

    length_type getCellSize() const
    {
        return impl_.cell_size();
    }

#if OBJECTMATRIX_USE_ITERATOR
    Enumerator<const Generators::result_type&>* iterneighbors(
            const sphere<length_type>& sphere)
    {
        return Generators::enumerate_neighbors(*this, sphere);
    }

    Enumerator<const Generators::result_type&>* iterneighbors_cyclic(
            const sphere<length_type>& sphere)
    {
        return Generators::enumerate_neighbors_cyclic(*this, sphere);
    }
#endif /* OBJECTMATRIX_USE_ITERATOR */

    typename impl_type::const_iterator __iter__begin() const
    {
        return impl_.begin();
    }

    typename impl_type::const_iterator __iter__end()
    {
        return impl_.end();
    }

    typename boost::transform_iterator<
        select1st<const typename impl_type::const_iterator::value_type>,
        typename impl_type::const_iterator> iterkeys_begin()
    {
        return boost::make_transform_iterator(impl_.begin(),
            select1st<const typename impl_type::const_iterator::value_type>());
    }

    typename boost::transform_iterator<
        select1st<const typename impl_type::const_iterator::value_type>,
        typename impl_type::const_iterator> iterkeys_end()
    {
        return boost::make_transform_iterator(impl_.end(),
            select1st<const typename impl_type::const_iterator::value_type>());
    }

    boost::shared_ptr<typename Builders::result_type>
    getNeighborsWithinRadiusNoSort(const position_type& pos, length_type radius)
    {
        typename Builders::distance_array_type::allocator_type alloc;

        if (radius >= impl_.cell_size() / 2)
        {
            throw std::runtime_error("Radius must be smaller than the half of the cell size");
        }

        boost::shared_ptr<typename Builders::result_type> retval(
            new typename Builders::result_type(
                typename boost::tuples::element<0, typename Builders::result_type>::type(),
                typename boost::tuples::element<1, typename Builders::result_type>::type(alloc)));
        Builders::build_neighbors_array_cyclic(*retval, *this,
                ::sphere<length_type>( pos, radius ) );

        // give away the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    boost::shared_ptr<typename Builders::result_type>
    getNeighborsCyclicNoSort(const position_type& pos)
    {
        typename Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<typename Builders::result_type> retval(
            new typename Builders::result_type(
                typename boost::tuples::element<0, typename Builders::result_type>::type(),
                typename boost::tuples::element<1, typename Builders::result_type>::type(alloc)));
        Builders::build_all_neighbors_array_cyclic(*retval, *this, pos);

        // give away the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    const bool contains(const key_type& k)
    {
        return impl_.find(k) != impl_.end();
    }

    const mapped_type& __getitem__(const key_type& k)
    {
        typename impl_type::iterator i(impl_.find(k));
        if (i == impl_.end())
        {
            PyErr_SetObject(PyExc_KeyError,
                    boost::python::incref(boost::python::object(k).ptr()));
            boost::python::throw_error_already_set();
        }

        return (*i).second;
    }

    bool __setitem__(const key_type& key, mapped_type const& val)
    {
        return impl_.insert(typename impl_type::value_type(key, val)).second;
    }

    void __delitem__(const key_type& key)
    {
        impl_.erase(key);
    }

    void check()
    {
        // do nothing
    }

    inline static void __register_class(char const* class_name)
    {
        using namespace boost::python;

        util::register_tuple_converter<typename impl_type::value_type>();

#if OBJECTMATRIX_USE_ITERATOR
        util::register_tuple_converter<typename Generators::result_type>();

        util::register_enumerator<
                const Generators::result_type&,
                return_value_policy<copy_const_reference> >(
                "ObjectContainer_NeighborIterator");
#endif /* OBJECTMATRIX_USE_ITERATOR */

        util::register_multi_array_converter<
            typename boost::tuples::element<0, typename Builders::result_type>::type>();
        util::register_multi_array_converter<
            typename boost::tuples::element<1, typename Builders::result_type>::type>();
        // the following conversion is the same as the previous
        // util::register_multi_array_converter<boost::tuples::element<2, Builders::result_type>::type>();

        util::register_tuple_converter<typename Builders::result_type>();

        class_<ObjectContainer>(class_name, init<length_type, size_type>())
#if OBJECTMATRIX_USE_ITERATOR
            .def("iterneighbors", &ObjectContainer::iterneighbors,
                    return_value_policy<manage_new_object>())
            .def("iterneighbors_cyclic", &ObjectContainer::iterneighbors_cyclic,
                    return_value_policy<manage_new_object>())
#endif /* OBJECTMATRIX_USE_ITERATOR */
            .add_property("cellSize", &ObjectContainer::getCellSize)
            .add_property("worldSize", &ObjectContainer::getWorldSize)
            .add_property("matrixSize", &ObjectContainer::getMatrixSize)
            .def("getNeighborsWithinRadiusNoSort", &ObjectContainer::getNeighborsWithinRadiusNoSort)
            .def("getNeighborsCyclicNoSort", &ObjectContainer::getNeighborsCyclicNoSort)
            .def("getNeighborsNoSort", &ObjectContainer::getNeighborsCyclicNoSort)
            .def("__len__", &ObjectContainer::__len__)
            .def("__iter__", range(
                    &ObjectContainer::__iter__begin,
                    &ObjectContainer::__iter__end))
            .def("iterkeys", range(
                    &ObjectContainer::iterkeys_begin,
                    &ObjectContainer::iterkeys_end))
            .def("contains", &ObjectContainer::contains)
            .def("__setitem__", &ObjectContainer::__setitem__)
            .def("__getitem__", &ObjectContainer::__getitem__,
                    return_value_policy<copy_const_reference>())
            .def("__delitem__", &ObjectContainer::__delitem__)
            .def("check", &ObjectContainer::check)
            ;
    }

private:
    impl_type impl_;
};

} // namespace peer

#endif /* OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP */
