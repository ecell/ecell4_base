#ifndef OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP
#define OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP

#include "../../../config.h"

#include <functional>
#include <string>
#include <vector>

#if HAVE_UNORDERED_MAP
#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#elif HAVE_EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <map>
#endif /* HAVE_UNORDERED_MAP */

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

#include "position.hpp"
#include "filters.hpp"
#include "sphere.hpp"
#include "object_container.hpp"

#if OBJECTMATRIX_USE_ITERATOR
#include <boost/coroutine/generator.hpp>
#include <boost/bind.hpp>
#include "peer/generator_support.hpp"
#endif /* OBJECTMATRIX_USE_ITERATOR */

namespace peer {

template<typename Tkey_, typename Tval_>
struct get_mapper_mf
{
#if HAVE_UNORDERED_MAP
    typedef std::unordered_map<Tkey_, Tval_> type;
#elif HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<Tkey_, Tval_> type;
#elif HAVE_EXT_HASH_MAP
    typedef __gnu_cxx::hash_map<Tkey_, Tval_> type;
#else 
    typedef std::map<Tkey_, Tval_> type;
#endif
};

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

class ObjectContainer
{
public:
    typedef boost::python::object key_type;
    typedef ::object_container<double, key_type, get_mapper_mf> impl_type;
    typedef impl_type::mapped_type mapped_type;
    typedef impl_type::position_type position_type;
    typedef impl_type::length_type length_type;
    typedef impl_type::size_type size_type;
    typedef impl_type::matrix_type::size_type matrix_size_type;

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
        typedef std::vector<impl_type::iterator> sphere_ref_array_type;
        typedef std::vector<impl_type::key_type> key_array_type;
        typedef std::vector<double, util::pyarray_backed_allocator<double> >
                distance_array_type;
        //typedef boost::tuple<sphere_ref_array_type, distance_array_type>
//                result_type;
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
                : //sa_(boost::get<0>(result)),
                ka_(boost::get<0>(result)),
                da_(boost::get<1>(result)) {}

            inline void operator()(impl_type::iterator i,
                    const position_type::value_type& d)
            {
                //sa_.push_back(i);
                ka_.push_back((*i).first);
                da_.push_back(d);
            }

        private:
            //sphere_ref_array_type& sa_;
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
                : //sa_(boost::get<0>(result)),
                ka_(boost::get<0>(result)),
                      da_(boost::get<1>(result)),
                      pos_(pos) {}

            inline void operator()(impl_type::iterator i)
            {
                //sa_.push_back(i);
                ka_.push_back((*i).first);
                da_.push_back(pos_.distance((*i).second.position) 
                              - (*i).second.radius);

            }

            inline void operator()(impl_type::iterator i,
                    const position_type& d)
            {
                //sa_.push_back(i);
                ka_.push_back((*i).first);
                da_.push_back(pos_.distance((*i).second.position + d)
                              - (*i).second.radius);
            }

        private:
            //sphere_ref_array_type& sa_;
            key_array_type& ka_;
            distance_array_type& da_;
            position_type pos_;
        };


    public:
        inline static void
        build_neighbors_array(result_type& retval,
                              ObjectContainer& cntnr, const mapped_type& sphere)
        {
            collector col(retval);
            take_neighbor(cntnr.impl_, col, sphere);
        }

        inline static void
        build_neighbors_array_cyclic(result_type& retval,
                ObjectContainer& cntnr, const mapped_type& sphere)
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
    ObjectContainer() {}

    ObjectContainer(length_type world_size, matrix_size_type size)
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

#if OBJECTMATRIX_USE_ITERATOR
    Enumerator<const Generators::result_type&>* iterneighbors(
            const sphere<double>& sphere)
    {
        return Generators::enumerate_neighbors(*this, sphere);
    }

    Enumerator<const Generators::result_type&>* iterneighbors_cyclic(
            const sphere<double>& sphere)
    {
        return Generators::enumerate_neighbors_cyclic(*this, sphere);
    }
#endif /* OBJECTMATRIX_USE_ITERATOR */

    boost::shared_ptr<Builders::result_type>
    //neighbors_array(const sphere& sphere)
    neighbors_array(const position_type& pos, const double radius)
    {
        Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<Builders::result_type> retval(
            new Builders::result_type(
                boost::tuples::element<0, Builders::result_type>::type(),
                boost::tuples::element<1, Builders::result_type>::type(alloc)));
        Builders::build_neighbors_array(*retval, *this,
                                        sphere<double>( pos, radius ) );

        // take over the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    boost::shared_ptr<Builders::result_type>
    //neighbors_array_cyclic(const sphere& sphere)
    neighbors_array_cyclic(const position_type& pos, const double radius)
    {
        Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<Builders::result_type> retval(
            new Builders::result_type(
                boost::tuples::element<0, Builders::result_type>::type(),
                boost::tuples::element<1, Builders::result_type>::type(alloc)));
        Builders::build_neighbors_array_cyclic(*retval, *this,
                                               sphere<double>( pos, radius ) );

        // take over the ownership of the arrays to the Numpy facility
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

        // take over the ownership of the arrays to the Numpy facility
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

        // take over the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    const bool contains( const key_type& k )
    {
        impl_type::iterator i(impl_.find(k));
        if (i == impl_.end())
        {
            return false;
        }
        return true;
    }

    const boost::tuple<position_type,double> get( const key_type& k )
    {
        impl_type::iterator i(impl_.find(k));
        if (i == impl_.end())
        {
            // FIXME: throw exception
            return boost::make_tuple(position_type(), 0.);
        }
        return boost::make_tuple( (*i).second.position, (*i).second.radius );
    }

    void insert( const key_type& key, const position_type& p, const double r )
    {
        impl_.insert(impl_type::value_type(key, sphere<double>(p,r)));
    }

    void update( const key_type& key, const position_type& p, const double r )
    {
        impl_.insert(impl_type::value_type(key, sphere<double>(p,r)));
    }


    void erase(const key_type& key)
    {
        impl_.erase(key);
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


#if OBJECTMATRIX_USE_ITERATOR
        util::register_tuple_converter<Generators::result_type>();

        util::register_enumerator<
                const Generators::result_type&,
                return_value_policy<copy_const_reference> >(
                "ObjectContainer_NeighborIterator");
#endif /* OBJECTMATRIX_USE_ITERATOR */

        util::register_multi_array_converter<
            boost::tuples::element<0, Builders::result_type>::type>();
        util::register_multi_array_converter<
            boost::tuples::element<1, Builders::result_type>::type>();
        // the following conversion is the same as the previous
        // util::register_multi_array_converter<boost::tuples::element<2, Builders::result_type>::type>();

        util::register_tuple_converter<Builders::result_type>();

        class_<ObjectContainer>("ObjectContainer")
            .def(init<length_type, matrix_size_type>())
#if OBJECTMATRIX_USE_ITERATOR
            .def("iterneighbors", &ObjectContainer::iterneighbors,
                    return_value_policy<manage_new_object>())
            .def("iterneighbors_cyclic", &ObjectContainer::iterneighbors_cyclic,
                    return_value_policy<manage_new_object>())
#endif /* OBJECTMATRIX_USE_ITERATOR */
            .add_property("cell_size", &ObjectContainer::cell_size)
            .add_property("world_size", &ObjectContainer::world_size)
            .add_property("matrix_size", &ObjectContainer::matrix_size)
            .def("neighbors_array", &ObjectContainer::neighbors_array)
            .def("neighbors_array_cyclic", &ObjectContainer::neighbors_array_cyclic)
            .def("all_neighbors_array", &ObjectContainer::all_neighbors_array)
            .def("all_neighbors_array_cyclic", &ObjectContainer::all_neighbors_array_cyclic)
            .def("size", &ObjectContainer::size)
            .def("contains", &ObjectContainer::contains)
            .def("insert", &ObjectContainer::insert)
            .def("update", &ObjectContainer::update)
            .def("get", &ObjectContainer::get)
            .def("erase", &ObjectContainer::erase);
    }

private:
    impl_type impl_;
};

} // namespace peer

#endif /* OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP */
