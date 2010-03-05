#ifndef OBJECTMATRIX_PEER_MATRIXSPACE_HPP
#define OBJECTMATRIX_PEER_MATRIXSPACE_HPP

#include <cstddef>
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
#include <boost/python/back_reference.hpp>
#include <boost/type_traits/alignment_of.hpp>
#include <numpy/arrayobject.h>

#include "peer/utils.hpp"
#include "peer/numpy/type_mappings.hpp"

#include "peer/tuple_converters.hpp"
#include "peer/numpy/pyarray_backed_allocator.hpp"
#include "peer/numpy/ndarray_converters.hpp"
#include "peer/STLIteratorWrapper.hpp"

#include "utils/pair.hpp"
#include "utils/range.hpp"
#include "filters.hpp"
#include "Shape.hpp"

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

    template<typename Tlength_>
    class CollectorResultConverter
    {
    public:
        typedef Tlength_ length_type;
        typedef std::pair<PyObject*, length_type> result_element;
        typedef std::vector<result_element, util::pyarray_backed_allocator<result_element> > result_type;

    private:
        struct to_ndarray_converter
        {
            static PyObject* convert(const result_type& val)
            {
                const npy_intp dims[1] = { val.size() };
                boost::python::incref(reinterpret_cast<PyObject*>(result_type_descr_));
                PyObject* retval = PyArray_NewFromDescr(&PyArray_Type,
                        result_type_descr_,
                        1, const_cast<npy_intp*>(dims), NULL,
                        &const_cast<result_type&>(val)[0],
                        NPY_CARRAY, NULL);
                if (!retval)
                    return NULL;
                reinterpret_cast<PyArrayObject*>(retval)->flags |= NPY_OWNDATA;
                return retval;
            }
        };


    public:
        static void __register_converter()
        {
            if (!result_type_descr_)
            {
                init_result_type_descr();
                boost::python::to_python_converter<result_type, to_ndarray_converter>();
            }
        }

    private:
        static void init_result_type_descr()
        {
            namespace py = boost::python;
            py::dict fields;
            py::str _pair("pair");
            py::str _distance("distance");
            fields[_pair] = py::make_tuple(
                py::object(
                    py::detail::new_reference(
                        reinterpret_cast<PyObject*>(PyArray_DescrFromType(
                        util::get_numpy_typecode<
                            typename result_element::first_type>::value)))),
                offsetof(result_element, first));
            fields[_distance] = py::make_tuple(
                py::object(
                    py::detail::new_reference(
                        reinterpret_cast<PyObject*>(PyArray_DescrFromType(
                            util::get_numpy_typecode<
                                typename result_element::second_type>::value)))),
                offsetof(result_element, second));
            result_type_descr_ = PyArray_DescrNewFromType(PyArray_VOID);
            result_type_descr_->hasobject = NPY_ITEM_HASOBJECT; 
            result_type_descr_->fields = py::incref(fields.ptr());
            result_type_descr_->names = py::incref(py::make_tuple(_pair, _distance).ptr());
            result_type_descr_->elsize = sizeof(result_element);
            result_type_descr_->alignment = boost::alignment_of<std::size_t>::value;
        }

    private:
        static PyArray_Descr* result_type_descr_;
    };
};

template<typename Tlength_>
PyArray_Descr* MatrixSpaceBase::CollectorResultConverter<Tlength_>::result_type_descr_(0);

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
        typedef CollectorResultConverter<length_type> collector_result_converter_type; 
        typedef typename collector_result_converter_type::result_element result_element;
        typedef typename collector_result_converter_type::result_type result_type;
        typedef typename remove_const_first<typename impl_type::value_type>::type value_type;

        struct collector
        {
        public:
            inline collector(typename Builders::result_type& result)
                : result_(result) {}

            inline void operator()(typename impl_type::iterator const& i,
                    const typename position_type::value_type& d) const
            {
                result_.push_back(result_element(
                    boost::python::incref(boost::python::object(*i).ptr()), d));
            }

            inline void operator()(typename impl_type::const_iterator const& i,
                    const typename position_type::value_type& d) const
            {
                result_.push_back(result_element(
                    boost::python::incref(boost::python::object(*i).ptr()), d));
            }

        private:
            typename Builders::result_type& result_;
        };

        template<typename TdistFun_>
        struct all_neighbors_collector
        {
        public:
            typedef TdistFun_ distance_functor;
        public:
            inline all_neighbors_collector(typename Builders::result_type& result,
                    const distance_functor& distance)
                : result_(result),
                  distance_(distance) {}

            inline void operator()(typename impl_type::iterator i) const
            {
                result_.push_back(result_element(
                    boost::python::incref(boost::python::object(*i).ptr()),
                    distance_(shape((*i).second))));
            }

            inline void operator()(typename impl_type::const_iterator const& i) const
            {
                result_.push_back(result_element(
                    boost::python::incref(boost::python::object(*i).ptr()),
                    distance_(shape((*i).second))));
            }

            inline void operator()(typename impl_type::iterator i,
                    const position_type& d) const
            {
                result_.push_back(result_element(
                    boost::python::incref(boost::python::object(*i).ptr()),
                    distance_(offset(shape((*i).second), d))));
            }

            inline void operator()(typename impl_type::const_iterator const& i,
                    const position_type& d) const
            {
                result_.push_back(result_element(
                    boost::python::incref(boost::python::object(*i).ptr()),
                    distance_(offset(shape((*i).second), d))));
            }

        private:
            typename Builders::result_type& result_;
            distance_functor const& distance_;
        };

        struct distance_comparator:
                public std::binary_function<result_element, result_element, bool> {

            bool operator()(result_element const& lhs, result_element const& rhs) const
            {
                return lhs.second < rhs.second;
            }
        };

        struct distance_calculator
        {
            distance_calculator(position_type const& pos): pos_(pos) {}

            template<typename Tshape_>
            length_type operator()(Tshape_ const& shape) const
            {
                return distance(shape, pos_);
            }

        private:
            position_type pos_;
        };

        struct cyclic_distance_calculator
        {
            cyclic_distance_calculator(position_type const& pos, length_type const& world_size): pos_(pos), world_size_(world_size) {}

            template<typename Tshape_>
            length_type operator()(Tshape_ const& shape) const
            {
                return distance_cyclic(shape, pos_, world_size_);
            }

        private:
            position_type pos_;
            length_type world_size_;
        };

    public:
        inline static void
        build_neighbors_array(result_type& retval,
                              impl_type const& cntnr,
                              const ::Sphere<length_type>& sphere)
        {
            collector col(retval);
            take_neighbor(cntnr, col, sphere);
            std::sort(retval.begin(), retval.end(), distance_comparator());
        }

        inline static void
        build_neighbors_array_cyclic(result_type& retval,
                impl_type const& cntnr, const ::Sphere<length_type>& sphere)
        {
            collector col(retval);
            take_neighbor_cyclic(cntnr, col, sphere);
            std::sort(retval.begin(), retval.end(), distance_comparator());
        }

        inline static void
        build_all_neighbors_array(result_type& retval,
                impl_type const& cntnr, const position_type& pos)
        {
            distance_calculator distance(pos);
            all_neighbors_collector<distance_calculator> col(retval, distance);
            cntnr.each_neighbor(cntnr.index(pos), col);
            std::sort(retval.begin(), retval.end(), distance_comparator());
        }

        inline static void
        build_all_neighbors_array_cyclic(result_type& retval,
                impl_type const& cntnr, const position_type& pos)
        {
            cyclic_distance_calculator distance(pos, cntnr.world_size());
            all_neighbors_collector<cyclic_distance_calculator> col(retval, distance);
            cntnr.each_neighbor_cyclic(cntnr.index(pos), col);
            std::sort(retval.begin(), retval.end(), distance_comparator());
        }

        static void __register_converter()
        {
            collector_result_converter_type::__register_converter();
            util::register_tuple_converter<value_type>();
            util::register_tuple_converter<typename impl_type::value_type>();
        }

    private:
        Builders() {}

    private:
    };

    typedef boost::transform_iterator<
        select1st<const typename impl_type::const_iterator::value_type>,
        typename impl_type::const_iterator> key_iterator;

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

    static boost::python::object __iter__(boost::python::back_reference<MatrixSpace const&> self)
    {
        using namespace boost::python;
        return object(handle<>(
            STLIteratorWrapper<typename impl_type::const_iterator,
                               boost::python::object>::create(
                std::make_pair(self.get().impl_.begin(),
                               self.get().impl_.end()),
                               self.source())));
    }

    static boost::python::object iterkeys(boost::python::back_reference<MatrixSpace const&> self)
    {
        using namespace boost::python;
        return object(handle<>(
            STLIteratorWrapper<key_iterator, boost::python::object>::create(
                make_transform_iterator_range(
                    std::make_pair(self.get().impl_.begin(),
                                   self.get().impl_.end()),
                    select1st<const typename impl_type::const_iterator::value_type>()),
                self.source())));
    }

    typename Builders::result_type
    get_neighbors_within_radius_cyclic(const position_type& pos,
                                       length_type radius)
    {
        typename Builders::result_type::allocator_type alloc;

        if (radius >= impl_.cell_size() / 2)
        {
            throw std::runtime_error("Radius must be smaller than the half of the cell size");
        }

        typename Builders::result_type retval(alloc);

        Builders::build_neighbors_array_cyclic(retval, impl_,
                ::Sphere<length_type>( pos, radius ) );

        // take over the ownership of the arrays to the Numpy facility
        alloc.giveup_ownership();
        return retval;
    }

    typename Builders::result_type
    get_neighbors_cyclic(const position_type& pos)
    {
        typename Builders::result_type::allocator_type alloc;
        typename Builders::result_type retval(alloc);

        Builders::build_all_neighbors_array_cyclic(retval, impl_, pos);

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

        Builders::__register_converter();

        util::register_tuple_converter<typename impl_type::value_type>();
        STLIteratorWrapper<typename impl_type::const_iterator,
                           boost::python::object>::__class_init__(
            "MatrixSpace.iterator");
        STLIteratorWrapper<key_iterator, boost::python::object>::__class_init__(
            "MatrixSpace.keyiterator");

        class_<MatrixSpace>(class_name, init<length_type, size_type>())
            .add_property("cell_size", &MatrixSpace::get_cell_size)
            .add_property("world_size", &MatrixSpace::get_world_size)
            .add_property("matrix_size", &MatrixSpace::get_matrix_size)
            .def("get_neighbors_within_radius_cyclic", &MatrixSpace::get_neighbors_within_radius_cyclic)
            .def("get_neighbors_within_radius", &MatrixSpace::get_neighbors_within_radius_cyclic)
            .def("get_neighbors_cyclic", &MatrixSpace::get_neighbors_cyclic)
            .def("get_neighbors", &MatrixSpace::get_neighbors_cyclic)
            .def("__len__", &MatrixSpace::__len__)
            .def("__iter__", &MatrixSpace::__iter__)
            .def("iterkeys", &MatrixSpace::iterkeys)
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
