#ifndef OBJECTMATRIX_PEER_RANGE_CONVERTERS_HPP
#define OBJECTMATRIX_PEER_RANGE_CONVERTERS_HPP

#include <utility>
#include <Python.h>
#include <listobject.h>
#include <tupleobject.h>
#include <abstract.h>
#include <boost/python.hpp>
#include <boost/optional.hpp>
#include <boost/range/size.hpp>
#include <boost/range/const_iterator.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/end.hpp>
#include <boost/range/value_type.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include "peer/utils.hpp"
#include "generator.hpp"

namespace peer {

namespace util
{
    namespace detail
    {
        template<typename Trange_, typename Tpolicy_>
        struct range_to_pyseq_converter
        {
            typedef Trange_ native_type;
            typedef Tpolicy_ policy_type;
            
            static PyObject* convert(const native_type& p)
            {
                using namespace boost::python;
                PyObject* retval = policy_type::create(boost::size(p));
                Py_ssize_t idx = 0;
                for (typename boost::range_const_iterator<native_type>::type i(boost::begin(p)), e(boost::end(p)); i != e; ++i, ++idx)
                {
                    policy_type::set(retval, idx, incref(object(*i).ptr()));
                }
                return retval;
            }
        };

        struct tuple_policy
        {
            static PyObject* create(Py_ssize_t size)
            {
                return PyTuple_New(size);
            }

            static void set(PyObject* tuple, Py_ssize_t idx, PyObject* obj)
            {
                PyTuple_SET_ITEM(tuple, idx, obj);
            }
        };

        struct list_policy
        {
            static PyObject* create(Py_ssize_t size)
            {
                return PyList_New(size);
            }

            static void set(PyObject* tuple, Py_ssize_t idx, PyObject* obj)
            {
                PyList_SET_ITEM(tuple, idx, obj);
            }
        };

        template<typename Trange_>
        struct range_to_pytuple_converter
            : public range_to_pyseq_converter<Trange_, tuple_policy>
        {
            typedef Trange_ native_type;
        };

        template<typename Trange_>
        struct range_to_pylist_converter
            : public range_to_pyseq_converter<Trange_, list_policy>
        {
            typedef Trange_ native_type;
        };

        template<typename Trange_>
        struct pyiterable_to_range_converter
        {
            typedef Trange_ native_type;
 
            static void* convertible(PyObject* pyo)
            {
                PyObject* const retval(PyObject_GetIter(pyo));
                if (!retval)
                {
                    PyErr_Clear();
                    return 0;
                }
                return retval;
            }

            static void construct(PyObject* pyo,
                                  boost::python::converter::rvalue_from_python_stage1_data* data)
            {
                void* storage(reinterpret_cast<
                    boost::python::converter::rvalue_from_python_storage<native_type>*>(data)->storage.bytes);
                boost::python::handle<> iter(
                        reinterpret_cast<PyObject*>(data->convertible));

                data->convertible = new (storage) native_type();
                native_type& retval(*reinterpret_cast<native_type*>(data->convertible));
                for (;;)
                {
                    boost::python::handle<> i(
                            boost::python::allow_null(PyIter_Next(iter.get())));
                    if (!i)
                    {
                        if (PyErr_Occurred())
                        {
                            boost::python::throw_error_already_set();
                        }
                        break;
                    }
                    retval.insert(boost::end(retval),
                        boost::python::extract<
                            typename boost::range_value<native_type>::type>(
                                i.get())());
                }
            }
        };

        template<typename Trange_, std::size_t N_>
        struct pyiterable_to_ra_container_converter
        {
            typedef Trange_ native_type;
 
            static void* convertible(PyObject* pyo)
            {
                PyObject* const retval(PyObject_GetIter(pyo));
                if (!retval)
                {
                    PyErr_Clear();
                    return 0;
                }
                return retval;
            }

            static void construct(PyObject* pyo,
                                  boost::python::converter::rvalue_from_python_stage1_data* data)
            {
                void* storage(reinterpret_cast<
                    boost::python::converter::rvalue_from_python_storage<native_type>*>(data)->storage.bytes);
                boost::python::handle<> iter(
                        reinterpret_cast<PyObject*>(data->convertible));

                data->convertible = new (storage) native_type();
                native_type& retval(*reinterpret_cast<native_type*>(data->convertible));
                std::size_t idx(0);
                for (;;)
                {
                    boost::python::handle<> i(boost::python::allow_null(PyIter_Next(iter.get())));
                    if (!i)
                    {
                        if (PyErr_Occurred())
                        {
                            boost::python::throw_error_already_set();
                        }
                        break;
                    }
                    if (idx >= N_)
                    {
                        PyErr_Format(PyExc_ValueError, "iterable generated more than %zd items", N_);
                        boost::python::throw_error_already_set();
                    }
                    retval[idx++] = boost::python::extract<
                            typename boost::range_value<native_type>::type>(
                                i.get())();
                }
                if (idx < N_)
                {
                    PyErr_Format(PyExc_ValueError, "iterable generated less than %zd items", N_);
                    boost::python::throw_error_already_set();
                }
            }
        };

        template<typename Tvalue_>
        class pyiterator_generator: public abstract_limited_generator<Tvalue_>
        {
        public:
            typedef Tvalue_ value_type;

        public:
            pyiterator_generator(boost::python::handle<> iter)
                : iter_(iter), advanced_(false) {}

            virtual ~pyiterator_generator() {}

            virtual bool valid() const
            {
                const_cast<pyiterator_generator*>(this)->fetch();
                return last_;
            }

            virtual value_type operator()()
            {
                fetch();
                if (!last_)
                {
                    return value_type();
                }
                advanced_ = false;
                return boost::python::extract<value_type>(last_.get())();
            }

            bool operator==(pyiterator_generator const& rhs) const
            {
                return (!last_ && !rhs.last_) || iter_ == rhs.iter_;
            }

        private:
            void fetch()
            {
                if (iter_ && !advanced_)
                {
                    last_ = boost::python::handle<>(
                            boost::python::allow_null(
                                PyIter_Next(iter_.get())));
                    if (!last_)
                    {
                        iter_.reset();
                    }
                    advanced_ = true;
                }
            }

        protected:
            boost::python::handle<> iter_; 
            bool advanced_;
            boost::python::handle<> last_;
        };

        template<typename Tvalue_>
        class pyseq_iterator
            : public boost::iterator_facade<
                pyseq_iterator<Tvalue_>, Tvalue_,
                boost::random_access_traversal_tag,
                Tvalue_,
                Py_ssize_t>
        {
        public:
            typedef Py_ssize_t difference_type;
            typedef Tvalue_ reference;

        public:
            pyseq_iterator(boost::python::object seq, difference_type idx = 0)
                : seq_(seq), idx_(idx) {}

            reference dereference() const
            {
                return boost::python::extract<Tvalue_>(
                    boost::python::handle<>(
                        PySequence_GetItem(seq_.ptr(), idx_)).get())();
            }

            bool equal(pyseq_iterator const& rhs) const
            {
                return seq_ == rhs.seq_ && idx_ == rhs.idx_;
            }

            void increment()
            {
                ++idx_;
            }

            void decrement()
            {
                --idx_;
            }

            void advance(difference_type n)
            {
                idx_ += n; 
            }

            difference_type distance_to(pyseq_iterator const& rhs) const
            {
                return rhs.idx_ - idx_;
            }

        protected:
            boost::python::object seq_;
            difference_type idx_;
        };
    }

    template<typename Tvalue_>
    struct py_range_wrapper
    {
        typedef detail::pyiterator_generator<Tvalue_> generator_type;

        struct generator_holder
        {
            bool operator==(generator_holder const& rhs) const
            {
                return impl_ == rhs.impl_;
            }

            operator bool() const
            {
                return impl_.is_initialized();
            }

            generator_type* operator->() const
            {
                return boost::get_pointer(impl_);
            }

            generator_type& operator*() const
            {
                return boost::get(impl_);
            }

            generator_holder(generator_type const& impl): impl_(impl) {}

            generator_holder(): impl_() {}

            mutable boost::optional<generator_type> impl_;
        };

        typedef generator_iterator<generator_type, generator_holder> iterator;
        typedef typename boost::remove_reference<
                typename generator_type::result_type> value_type;
        typedef typename generator_type::result_type reference;
        typedef iterator const_iterator;

        py_range_wrapper(boost::python::object obj): obj_(obj) {}

        std::size_t size() const
        {
            return PyObject_Size(obj_.ptr());
        }

        iterator begin() const
        {
            return iterator(generator_holder(generator_type(
                boost::python::handle<>(PyObject_GetIter(obj_.ptr())))));
        }

        iterator end() const
        {
            return iterator();
        }

    protected:
        boost::python::object obj_;
    };

    template<typename Tvalue_>
    struct py_range_wrapper_converter
    {
        typedef py_range_wrapper<Tvalue_> native_type;

        static void* convertible(PyObject* pyo)
        {
            if (!PyType_HasFeature(Py_TYPE(pyo), Py_TPFLAGS_HAVE_ITER))
            {
                return 0;
            }
            return pyo;
        }

        static void construct(PyObject* pyo,
                              boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            using namespace boost::python;

            void* storage(reinterpret_cast<
                converter::rvalue_from_python_storage<native_type>*>(data)->storage.bytes);
            data->convertible = new (storage) native_type(
                boost::python::object(boost::python::borrowed(pyo)));
        }
    };

    template<typename Tvalue_>
    struct pyseq_range_wrapper_converter
    {
        typedef boost::iterator_range<detail::pyseq_iterator<Tvalue_> > native_type;

        static void* convertible(PyObject* pyo)
        {
            if (!PySequence_Check(pyo))
            {
                return 0;
            }
            return pyo;
        }

        static void construct(PyObject* pyo,
                              boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            using namespace boost::python;

            void* storage(reinterpret_cast<
                converter::rvalue_from_python_storage<native_type>*>(data)->storage.bytes);
            data->convertible = new (storage) native_type(
                object(boost::python::detail::new_reference(pyo)));
        }
    };

} // namespace util

} // namespace peer

namespace boost
{

template<typename Tvalue>
inline typename boost::range_difference<peer::util::py_range_wrapper<Tvalue> >::type size(peer::util::py_range_wrapper<Tvalue> const& w)
{
    return w.size();
}

} // namespace boost

namespace peer {

namespace util
{

    template<typename Trange_>
    inline void register_range_to_tuple_converter()
    {
        static bool registered = false;
        if (!registered)
        {
            boost::python::to_python_converter<
                Trange_, detail::range_to_pytuple_converter<Trange_> >();
            registered = true;
        }
    }

    template<typename Trange_>
    inline void register_range_to_list_converter()
    {
        static bool registered = false;
        if (!registered)
        {
            boost::python::to_python_converter<
                Trange_, detail::range_to_pylist_converter<Trange_> >();
            registered = true;
        }
    }

    template<typename Trange_>
    inline void register_iterable_to_range_converter()
    {
        static bool registered = false;
        if (!registered)
        {
            to_native_converter<Trange_, detail::pyiterable_to_range_converter<Trange_> >();
            registered = true;
        }
    }

    template<typename Trange_, std::size_t N_>
    inline void register_iterable_to_ra_container_converter()
    {
        static bool registered = false;
        if (!registered)
        {
            to_native_converter<Trange_, detail::pyiterable_to_ra_container_converter<Trange_, N_> >();
            registered = true;
        }
    }

    template<typename Tvalue_>
    inline void register_py_range_wrapper_converter()
    {
        static bool registered = false;
        if (!registered)
        {
            to_native_converter<py_range_wrapper<Tvalue_>,
                                py_range_wrapper_converter<Tvalue_> >();
            registered = true;
        }
    }

} // namespace util

} // namespace peer

template<typename Tvalue>
inline bool valid(peer::util::detail::pyiterator_generator<Tvalue> const& gen)
{
    return gen.valid();
}

#endif /* OBJECTMATRIX_PEER_RANGE_CONVERTERS_HPP */
