#ifndef OBJECTMATRIX_PEER_MATRIXSPACE_HPP
#define OBJECTMATRIX_PEER_MATRIXSPACE_HPP

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

#include "filters.hpp"
#include "utils/pair.hpp"

namespace peer {

struct MatrixSpaceBase
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
class MatrixSpace: public MatrixSpaceBase
{
public:
    typedef Timpl_ impl_type;
    typedef typename impl_type::key_type key_type;
    typedef typename impl_type::mapped_type mapped_type;
    typedef typename impl_type::position_type position_type;
    typedef typename impl_type::length_type length_type;
    typedef typename impl_type::size_type size_type;
    typedef typename impl_type::matrix_type::size_type matrix_size_type;

    class Builders
    {
    public:
        typedef typename remove_const_first<typename impl_type::value_type>::type value_type;
        typedef boost::python::list value_array_type;
        typedef std::vector<length_type, util::pyarray_backed_allocator<length_type> >
                distance_array_type;
        typedef boost::tuple<value_array_type, distance_array_type>
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
                : sa_(boost::get<0>(result)),
                  da_(boost::get<1>(result)) {}

            inline void operator()(typename impl_type::iterator i,
                    const typename position_type::value_type& d)
            {
                sa_.append(boost::python::object(*i));
                da_.push_back(d);
            }

        private:
            value_array_type& sa_;
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
                : sa_(boost::get<0>(result)),
                  da_(boost::get<1>(result)),
                  pos_(pos) {}

            inline void operator()(typename impl_type::iterator i)
            {
                sa_.append(boost::python::object(*i));
                da_.push_back(distance(pos_, (*i).second.position()) 
                              - (*i).second.radius());

            }

            inline void operator()(typename impl_type::iterator i,
                    const position_type& d)
            {
                sa_.append(boost::python::object(*i));
                da_.push_back(distance(pos_, (*i).second.position() + d)
                              - (*i).second.radius());
            }

        private:
            value_array_type& sa_;
            distance_array_type& da_;
            position_type pos_;
        };


    public:
        inline static void
        build_neighbors_array(result_type& retval,
                              MatrixSpace& cntnr,
                              const ::Sphere<length_type>& sphere)
        {
            collector col(retval);
            take_neighbor(cntnr.impl_, col, sphere);
        }

        inline static void
        build_neighbors_array_cyclic(result_type& retval,
                MatrixSpace& cntnr, const ::Sphere<length_type>& sphere)
        {
            collector col(retval);
            take_neighbor_cyclic(cntnr.impl_, col, sphere);
        }

        inline static void
        build_all_neighbors_array(result_type& retval,
                MatrixSpace& cntnr, const position_type& pos)
        {
            all_neighbors_collector col(retval, pos);
            cntnr.impl_.each_neighbor(cntnr.impl_.index(pos), col);
        }

        inline static void
        build_all_neighbors_array_cyclic(result_type& retval,
                MatrixSpace& cntnr, const position_type& pos)
        {
            all_neighbors_collector col(retval, pos);
            cntnr.impl_.each_neighbor_cyclic(cntnr.impl_.index(pos), col);
        }

    private:
        Builders() {}
    };

public:
    MatrixSpace(typename impl_type::length_type world_size,
            typename impl_type::matrix_type::size_type size)
            : impl_(world_size, size) {}

    size_type __len__() const
    {
        return impl_.size();
    }

    size_type get_matrix_size() const
    {
        return impl_.matrix_size();
    }

    length_type get_world_size() const
    {
        return impl_.world_size();
    }

    length_type get_cell_size() const
    {
        return impl_.cell_size();
    }

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

    typename Builders::result_type
    get_neighbors_within_radius(const position_type& pos, length_type radius)
    {
        typename Builders::distance_array_type::allocator_type alloc;

        if (radius >= impl_.cell_size() / 2)
        {
            throw std::runtime_error("Radius must be smaller than the half of the cell size");
        }

        typename Builders::result_type retval(
            (typename boost::tuples::element<0, typename Builders::result_type>::type)typename boost::tuples::element<0, typename Builders::result_type>::type(),
            typename boost::tuples::element<1, typename Builders::result_type>::type(alloc));
        Builders::build_neighbors_array_cyclic(retval, *this,
                ::Sphere<length_type>( pos, radius ) );

        // take over the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    typename Builders::result_type
    get_neighbors_cyclic(const position_type& pos)
    {
        typename Builders::distance_array_type::allocator_type alloc;

       typename Builders::result_type retval(
            (typename boost::tuples::element<0, typename Builders::result_type>::type)typename boost::tuples::element<0, typename Builders::result_type>::type(),
            typename boost::tuples::element<1, typename Builders::result_type>::type(alloc));
        Builders::build_all_neighbors_array_cyclic(retval, *this, pos);

        // take over the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    const bool contains(const key_type& k)
    {
        return impl_.find(k) != impl_.end();
    }

    bool update(typename impl_type::value_type const& pair)
    {
        return impl_.update(pair).second;
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
        return impl_.update(typename impl_type::value_type(key, val)).second;
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

        util::register_tuple_converter<typename Builders::value_type>();

        //util::register_multi_array_converter<
        //    typename boost::tuples::element<0, typename Builders::result_type>::type>();
        util::register_multi_array_converter<
            typename boost::tuples::element<1, typename Builders::result_type>::type>();
        // the following conversion is the same as the previous
        // util::register_multi_array_converter<boost::tuples::element<2, Builders::result_type>::type>();

        util::register_tuple_converter<typename Builders::result_type>();

        class_<MatrixSpace>(class_name, init<length_type, size_type>())
            .add_property("cell_size", &MatrixSpace::get_cell_size)
            .add_property("world_size", &MatrixSpace::get_world_size)
            .add_property("matrix_size", &MatrixSpace::get_matrix_size)
            .def("get_neighbors_within_radius", &MatrixSpace::get_neighbors_within_radius)
            .def("get_neighbors_cyclic", &MatrixSpace::get_neighbors_cyclic)
            .def("get_neighbors", &MatrixSpace::get_neighbors_cyclic)
            .def("__len__", &MatrixSpace::__len__)
            .def("__iter__", range(
                    &MatrixSpace::__iter__begin,
                    &MatrixSpace::__iter__end))
            .def("iterkeys", range(
                    &MatrixSpace::iterkeys_begin,
                    &MatrixSpace::iterkeys_end))
            .def("contains", &MatrixSpace::contains)
            .def("update", &MatrixSpace::update)
            .def("__setitem__", &MatrixSpace::__setitem__)
            .def("__getitem__", &MatrixSpace::__getitem__,
                    return_value_policy<copy_const_reference>())
            .def("__delitem__", &MatrixSpace::__delitem__)
            .def("check", &MatrixSpace::check)
            ;
    }

private:
    impl_type impl_;
};

} // namespace peer

#endif /* OBJECTMATRIX_PEER_MATRIXSPACE_HPP */
