#ifndef OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP
#define OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP

#include <functional>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits.hpp>
#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <numpy/arrayobject.h>

#include "peer/utils.hpp"
#include "peer/numpy/type_mappings.hpp"

#include "peer/tuple_converters.hpp"
#include "peer/numpy/pyarray_backed_allocator.hpp"
#include "peer/numpy/ndarray_converters.hpp"
#include "peer/py_hash_support.hpp"

#include "filters.hpp"
#include "Sphere.hpp"
#include "../MatrixSpace.hpp"
#include "utils/get_mapper_mf.hpp"

namespace peer {

class MatrixSpace
{
public:
    typedef boost::python::object key_type;
    typedef ::MatrixSpace< Sphere<double>, key_type, get_mapper_mf> impl_type;
    typedef impl_type::mapped_type mapped_type;
    typedef impl_type::position_type position_type;
    typedef impl_type::length_type length_type;
    typedef impl_type::size_type size_type;
    typedef impl_type::matrix_type::size_type matrix_size_type;

    class Builders
    {
    public:
        typedef std::vector<impl_type::iterator> sphere_ref_array_type;
        typedef std::vector<impl_type::key_type> key_array_type;
        typedef std::vector<double, util::pyarray_backed_allocator<double> >
                distance_array_type;
        typedef boost::tuple<key_array_type, distance_array_type>
                result_type;

        struct collector: public std::binary_function<
                impl_type::reference, position_type::value_type, void>
        {
            typedef impl_type::iterator first_argument_type;
            typedef position_type::value_type second_argument_type;
            typedef void result_type;
        public:
            inline collector(Builders::result_type& result)
                : ka_(boost::get<0>(result)),
                  da_(boost::get<1>(result)) {}

            inline void operator()(impl_type::iterator i,
                    const position_type::value_type& d)
            {
                ka_.push_back((*i).first);
                da_.push_back(d);
            }

        private:
            key_array_type& ka_;
            distance_array_type& da_;
        };

        struct all_neighbors_collector: public std::binary_function<
                impl_type::reference, position_type::value_type, void>
        {
            typedef impl_type::iterator first_argument_type;
            typedef position_type::value_type second_argument_type;
            typedef void result_type;
        public:
            inline all_neighbors_collector(Builders::result_type& result,
                    const position_type& pos)
                : ka_(boost::get<0>(result)),
                  da_(boost::get<1>(result)),
                  pos_(pos) {}

            inline void operator()(impl_type::iterator i)
            {
                ka_.push_back((*i).first);
                da_.push_back(distance(pos_, (*i).second.position())
                              - (*i).second.radius());

            }

            inline void operator()(impl_type::iterator i,
                    const position_type& d)
            {
                ka_.push_back((*i).first);
                da_.push_back(distance(pos_, (*i).second.position() + d)
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
                              MatrixSpace& cntnr, const mapped_type& sphere)
        {
            collector col(retval);
            take_neighbor(cntnr.impl_, col, sphere);
        }

        inline static void
        build_neighbors_array_cyclic(result_type& retval,
                MatrixSpace& cntnr, const mapped_type& sphere)
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
    MatrixSpace() {}

    MatrixSpace(length_type world_size, matrix_size_type size)
        : impl_(world_size, size) {}

    size_type size() const
    {
        return impl_.size();
    }

    size_type matrix_size() const
    {
        return impl_.matrix_size();
    }

    length_type world_size() const
    {
        return impl_.world_size();
    }

    length_type cell_size() const
    {
        return impl_.cell_size();
    }

    boost::shared_ptr<Builders::result_type>
    neighbors_array(const position_type& pos, const double radius)
    {
        Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<Builders::result_type> retval(
            new Builders::result_type(
                boost::tuples::element<0, Builders::result_type>::type(),
                boost::tuples::element<1, Builders::result_type>::type(alloc)));
        Builders::build_neighbors_array(*retval, *this,
                                        Sphere<double>( pos, radius ) );

        // give away the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    boost::shared_ptr<Builders::result_type>
    neighbors_array_cyclic(const position_type& pos, const double radius)
    {
        Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<Builders::result_type> retval(
            new Builders::result_type(
                boost::tuples::element<0, Builders::result_type>::type(),
                boost::tuples::element<1, Builders::result_type>::type(alloc)));
        Builders::build_neighbors_array_cyclic(*retval, *this,
                                               Sphere<double>( pos, radius ) );

        // give away the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    boost::shared_ptr<Builders::result_type>
    all_neighbors_array(const position_type& pos)
    {
        Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<Builders::result_type> retval(
            new Builders::result_type(
                boost::tuples::element<0, Builders::result_type>::type(),
                boost::tuples::element<1, Builders::result_type>::type(alloc)));
        Builders::build_all_neighbors_array(*retval, *this, pos);

        // give away the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    boost::shared_ptr<Builders::result_type>
    all_neighbors_array_cyclic(const position_type& pos)
    {
        Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<Builders::result_type> retval(
            new Builders::result_type(
                boost::tuples::element<0, Builders::result_type>::type(),
                boost::tuples::element<1, Builders::result_type>::type(alloc)));
        Builders::build_all_neighbors_array_cyclic(*retval, *this, pos);

        // give away the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    bool contains( const key_type& k )
    {
        impl_type::iterator i(impl_.find(k));
        if (i == impl_.end())
        {
            return false;
        }
        return true;
    }


    mapped_type const& get( const key_type& k )
    {
        impl_type::iterator i(impl_.find(k));
        if (i == impl_.end())
        {
            throw std::runtime_error( "key not found." );
        }

        return (*i).second;
    }


    void update( const key_type& key, const position_type& p, const length_type& r )
    {
        impl_.update(impl_type::value_type(key, Sphere<length_type>(p,r)));
    }


    bool erase(const key_type& key)
    {
        return impl_.erase(key);
    }

    impl_type::const_iterator __iter__begin() const
    {
        return impl_.begin();
    }

    impl_type::const_iterator __iter__end() const
    {
        return impl_.end();
    }

    operator impl_type&()
    {
        return impl_;
    }

    operator const impl_type&() const
    {
        return impl_;
    }

    inline static void __register_class()
    {
        using namespace boost::python;

        util::register_tuple_converter<impl_type::value_type>();
        util::register_tuple_converter<boost::tuple<position_type, length_type> >();

#if OBJECTMATRIX_USE_ITERATOR
        util::register_tuple_converter<Generators::result_type>();

        util::register_enumerator<
                const Generators::result_type&,
                return_value_policy<copy_const_reference> >(
                "MatrixSpace_NeighborIterator");
#endif /* OBJECTMATRIX_USE_ITERATOR */

        util::register_multi_array_converter<
            boost::tuples::element<0, Builders::result_type>::type>();
        util::register_multi_array_converter<
            boost::tuples::element<1, Builders::result_type>::type>();
        // the following conversion is the same as the previous
        // util::register_multi_array_converter<boost::tuples::element<2, Builders::result_type>::type>();

        util::register_tuple_converter<Builders::result_type>();

        class_<MatrixSpace>("MatrixSpace")
            .def(init<length_type, matrix_size_type>())
#if OBJECTMATRIX_USE_ITERATOR
            .def("iterneighbors", &MatrixSpace::iterneighbors,
                    return_value_policy<manage_new_object>())
            .def("iterneighbors_cyclic", &MatrixSpace::iterneighbors_cyclic,
                    return_value_policy<manage_new_object>())
#endif /* OBJECTMATRIX_USE_ITERATOR */
            .add_property("cell_size", &MatrixSpace::cell_size)
            .add_property("world_size", &MatrixSpace::world_size)
            .add_property("matrix_size", &MatrixSpace::matrix_size)
            .def("neighbors_array", &MatrixSpace::neighbors_array)
            .def("neighbors_array_cyclic", &MatrixSpace::neighbors_array_cyclic)
            .def("all_neighbors_array", &MatrixSpace::all_neighbors_array)
            .def("all_neighbors_array_cyclic", &MatrixSpace::all_neighbors_array_cyclic)
            .def("size", &MatrixSpace::size)
            .def("contains", &MatrixSpace::contains)
            .def("update", &MatrixSpace::update)
            .def("get", &MatrixSpace::get,
                    return_value_policy<copy_const_reference>())
            .def("__iter__", range(
                    &MatrixSpace::__iter__begin,
                    &MatrixSpace::__iter__end))
            .def("erase", &MatrixSpace::erase);
    }

private:
    impl_type impl_;
};

} // namespace peer

#endif /* OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP */
