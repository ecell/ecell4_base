#ifndef OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP
#define OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP

#include <functional>
#include <string>
#include <vector>
#include <tr1/unordered_map>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits.hpp>
#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include "peer/tuple_converters.hpp"
#include "peer/numpy/pyarray_backed_allocator.hpp"
#include "peer/numpy/ndarray_converters.hpp"

#include "filters.hpp"
#include "object_container.hpp"

#if OBJECTMATRIX_USE_ITERATOR
#include <boost/coroutine/generator.hpp>
#include <boost/bind.hpp>
#include "peer/generator_support.hpp"
#endif /* OBJECTMATRIX_USE_ITERATOR */

#include "peer/Sphere.hpp"

namespace peer {

template<typename Tkey_, typename Tval_>
struct get_mapper_mf
{
    typedef std::tr1::unordered_map<Tkey_, Tval_> type;
};

template<typename Tval_>
struct get_mapper_mf<boost::python::object, Tval_>
{
    struct hasher: public std::unary_function<boost::python::object, std::size_t>
    {
        typedef boost::python::object argument_type;
        typedef std::size_t result_type;

        result_type operator()(const argument_type& arg) const
        {
            return static_cast<result_type>((long)PyObject_Hash(arg.ptr()));
        }
    };

    typedef std::tr1::unordered_map<boost::python::object, Tval_, hasher> type;
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
        typedef std::vector<double, util::pyarray_backed_allocator<double> >
                distance_array_type;
        typedef boost::tuple<sphere_ref_array_type, distance_array_type>
                result_type;

        struct collector: public std::binary_function<
                impl_type::reference, position_type::value_type, void>
        {
            typedef impl_type::iterator first_argument_type;
            typedef position_type::value_type second_argument_type;
            typedef void result_type;
        public:
            inline collector(Builders::result_type& result)
                    : sa_(boost::get<0>(result)),
                      da_(boost::get<1>(result)) {}

            inline void operator()(impl_type::iterator i,
                    const position_type::value_type& d)
            {
                sa_.push_back(i);
                da_.push_back(d);
            }

        private:
            sphere_ref_array_type& sa_;
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
                    : sa_(boost::get<0>(result)),
                      da_(boost::get<1>(result)),
                      pos_(pos) {}

            inline void operator()(impl_type::iterator i)
            {
                sa_.push_back(i);
                da_.push_back(pos_.distance((*i).second.position) 
                              - (*i).second.radius);

            }

            inline void operator()(impl_type::iterator i,
                    const position_type& d)
            {
                sa_.push_back(i);
                da_.push_back(pos_.distance((*i).second.position + d)
                              - (*i).second.radius);
            }

        private:
            sphere_ref_array_type& sa_;
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

    class SphereRef
    {
    public:
        typedef impl_type::mapped_type::value_type value_type;
        typedef impl_type::iterator impl_type;

    public:
        SphereRef(impl_type impl): impl_(impl) {}

        PyObject* __repr__() const
        {
            return util::pystr_from_repr(&(*impl_).second);
        }

        value_type _get_x() const
        {
            return (*impl_).second.position.x();
        }

        value_type _get_y() const
        {
            return (*impl_).second.position.y();
        }

        value_type _get_z() const
        {
            return (*impl_).second.position.z();
        }

        value_type _get_radius() const
        {
            return (*impl_).second.radius;
        }

        void _set_radius(value_type val)
        {
            (*impl_).second.radius = val;
        }

        boost::python::object _get_id() const
        {
            return (*impl_).first;
        }
   
        inline static void __register_class()
        {
            using namespace boost::python;

            class_<SphereRef>("SphereRef", no_init)
                .def("__repr__", &SphereRef::__repr__)
                .add_property("x", &SphereRef::_get_x)
                .add_property("y", &SphereRef::_get_y)
                .add_property("z", &SphereRef::_get_z)
                .add_property("radius",
                        &SphereRef::_get_radius,
                        &SphereRef::_set_radius)
                .add_property("id", &SphereRef::_get_id);
        }
 
    private:
        impl_type impl_;
    };

    struct iter_to_sphere_ref_converter
    {
        static PyObject* convert(const impl_type::iterator& i)
        {
            return boost::python::incref(boost::python::object(new SphereRef(i)).ptr());
        }
    };

public:
    ObjectContainer() {}

    ObjectContainer(length_type world_size, matrix_size_type size)
        : impl_(world_size, size) {}

    size_type __len__() const
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
            const Sphere& sphere)
    {
        return Generators::enumerate_neighbors(*this, sphere);
    }

    Enumerator<const Generators::result_type&>* iterneighbors_cyclic(
            const Sphere& sphere)
    {
        return Generators::enumerate_neighbors_cyclic(*this, sphere);
    }
#endif /* OBJECTMATRIX_USE_ITERATOR */

    boost::shared_ptr<Builders::result_type>
    neighbors_array(const Sphere& sphere)
    {
        Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<Builders::result_type> retval(
            new Builders::result_type(
                boost::tuples::element<0, Builders::result_type>::type(),
                boost::tuples::element<1, Builders::result_type>::type(alloc)));
        Builders::build_neighbors_array(*retval, *this,
                static_cast<const Sphere::impl_type&>(sphere));

        // take over the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    boost::shared_ptr<Builders::result_type>
    neighbors_array_cyclic(const Sphere& sphere)
    {
        Builders::distance_array_type::allocator_type alloc;

        boost::shared_ptr<Builders::result_type> retval(
            new Builders::result_type(
                boost::tuples::element<0, Builders::result_type>::type(),
                boost::tuples::element<1, Builders::result_type>::type(alloc)));
        Builders::build_neighbors_array_cyclic(*retval, *this,
                static_cast<const Sphere::impl_type&>(sphere));

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

    const bool __contains__(key_type k)
    {
        impl_type::iterator i(impl_.find(k));
        if (i == impl_.end())
        {
            return false;
        }
        return true;
    }

    SphereRef* __getitem__(key_type k)
    {
        impl_type::iterator i(impl_.find(k));
        if (i == impl_.end())
        {
            return NULL;
        }
        return new SphereRef(i);
    }

    void __setitem__(key_type key, Sphere* item)
    {
        impl_.insert(impl_type::value_type(key, *item));
    }

    void __delitem__(key_type key)
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

        SphereRef::__register_class();

        to_python_converter<impl_type::iterator,
                iter_to_sphere_ref_converter>();

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
            .def("__len__", &ObjectContainer::__len__)
            .def("__contains__", &ObjectContainer::__contains__)
            .def("__setitem__", &ObjectContainer::__setitem__)
            .def("__getitem__", &ObjectContainer::__getitem__,
                        return_value_policy<manage_new_object>())
            .def("__delitem__", &ObjectContainer::__delitem__);
    }

private:
    impl_type impl_;
};

} // namespace peer

#endif /* OBJECTMATRIX_PEER_OBJECTCONTAINER_HPP */
