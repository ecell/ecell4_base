#include <functional>
#include <stdexcept>
#include <string>
#include <complex>
#include <vector>
#include <limits>
#include <tr1/unordered_map>
#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/type_traits.hpp>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/object.hpp>
//#include <boost/coroutine/generator.hpp>
#include <Python.h>
#include <numpy/arrayobject.h>
#include "sphere.hpp"
#include "filters.hpp"
#include "object_container.hpp"

namespace peer {

template<typename T_>
struct Enumerator
{
    typedef T_ result_type;

    virtual ~Enumerator() {}

    virtual T_ next() = 0;
};

class StopIteration
{
public:
    StopIteration() {}
};


namespace util
{
    inline boost::python::object pass_through(boost::python::object o)
    {
        return o;
    }

    template<typename T_>
    struct EnumeratorWrapper
            : public Enumerator<T_>, boost::python::wrapper<Enumerator<T_> >
    {
        EnumeratorWrapper(const Enumerator<T_>& that)
            : Enumerator<T_>(that) {}

        T_ next()
        {
            return this->get_override("next")();
        }
    };

/*
    template<typename Tgen_>
    struct GeneratorEnumeratorAdapter
        : public Enumerator<typename boost::remove_pointer<
                typename Tgen_::result_type>::type&>
    {
        typedef typename boost::remove_pointer<
                typename Tgen_::result_type>::type& result_type;

        GeneratorEnumeratorAdapter(const Tgen_& gen): gen_(gen) {}

        virtual result_type next()
        {
            if (!*gen_)
                throw StopIteration();
            return *gen_();
        }

    private:
        Tgen_ gen_;
    };

    template<typename T_, typename Tpol_>
    void register_enumerator(const char* python_name)
    {
        boost::python::class_<EnumeratorWrapper<T_> >(python_name,
            boost::python::no_init)
            .def("next", &Enumerator<T_>::next, Tpol_())
            .def("__iter__", &pass_through);
    }
*/
    class pyarray_backed_allocator_base
    {
    public:
        typedef std::size_t size_type;
        typedef std::ptrdiff_t difference_type;

    protected:
        struct allocator_state
        {
            bool giveup_ownership;

        protected:
            std::size_t refcount;

        public:
            allocator_state(bool _giveup_ownership = false)
                : giveup_ownership(_giveup_ownership), refcount(1) {} 

            bool release()
            {
                return 0 == --refcount;
            }

            allocator_state& add_ref()
            {
                ++refcount;
                return *this;
            }
        };

    public:
        pyarray_backed_allocator_base(bool giveup_ownership)
            : state_(new allocator_state(giveup_ownership)) {}

        pyarray_backed_allocator_base(const pyarray_backed_allocator_base& that)
            : state_(&that.state_->add_ref()) {}

        ~pyarray_backed_allocator_base()
        {
            if (state_->release())
            {
                delete state_;
            }
        }

        void giveup_ownership()
        {
            
            state_->giveup_ownership = true;
        }

    protected:
        void* _nalloc(const size_type sz, const size_type n) const
        {
            if (static_cast<size_type>(static_cast<double>(sz) * n) !=
                    sz * n)
            {
                throw std::bad_alloc();
            }

            void* retval = PyArray_malloc(sz * n);
            if (!retval)
            {
                throw std::bad_alloc();
            }
            return retval;
        }

        // it's possible that "free" is previously defined as a
        // preprocessor macro.
        void _free(void* ptr) const
        {   
            if (ptr && !state_->giveup_ownership)
            {
                PyArray_free(ptr);
            }
        }

    protected:
        allocator_state* state_;
    };

    template<typename T_>
    class pyarray_backed_allocator: public pyarray_backed_allocator_base
    {
    public:
        typedef T_ value_type;
        typedef T_* pointer;
        typedef const T_* const_pointer;
        typedef T_& reference;
        typedef const T_& const_reference;

        template<typename Tother_>
        struct rebind
        {
            typedef pyarray_backed_allocator<Tother_> other;
        };

    public:
        pyarray_backed_allocator(bool giveup_ownership = false)
            : pyarray_backed_allocator_base(giveup_ownership) {}

        pyarray_backed_allocator(const pyarray_backed_allocator_base& that)
            : pyarray_backed_allocator_base(that) {}

        pointer address(reference r) const
        {
            return &r;
        }

        const_pointer address(const_reference r) const
        {
            return &r;
        }

        value_type* allocate(size_type n, const void* hint = 0) const
        {
            return reinterpret_cast<value_type*>(this->_nalloc(sizeof(T_), n));
        }

        void construct(T_* p, const T_& src) const
        {
            new(p) T_(src);
        }

        void destroy(T_* p) const
        {
            p->~T_(); // XXX: does this work for PODs?
        }

        size_type max_size() const
        {
            return std::numeric_limits<size_type>::max() / sizeof(T_);
        }

        void deallocate(T_* p, size_type n) const
        {
            this->_free(p);
        }

        bool operator==(const pyarray_backed_allocator_base& rhs)
        {
            return state_->giveup_ownership == rhs.state_->giveup_ownership;
        }

        bool operator!=(const pyarray_backed_allocator_base& rhs)
        {
            return !operator==(rhs);
        }
    };

    template<typename T_>
    inline PyObject* pystr_from_repr(const T_* obj)
    {
        std::ostringstream s;
        s << (*obj);
        std::string res(s.str());
        return PyString_FromStringAndSize(res.data(), res.size());
    }

    namespace detail
    {
        inline void _handle_stop_iteration_exc(const StopIteration& exc)
        {
            PyErr_SetNone(PyExc_StopIteration);
        }

        template<typename Ttcell_>
        inline void build_pytuple_from_tuple(PyObject* pyt, const Ttcell_& cell,
                Py_ssize_t idx = 0)
        {
            PyTuple_SetItem(pyt, idx,
                boost::python::incref(
                    boost::python::object(cell.get_head()).ptr()));
            build_pytuple_from_tuple(pyt, cell.get_tail(), idx + 1);
        }

        template<>
        inline void build_pytuple_from_tuple<boost::tuples::null_type>(
                PyObject*, const boost::tuples::null_type&, Py_ssize_t) {}

        template<typename Ttuple_>
        struct tuple_to_tuple_converter
        {
            typedef Ttuple_ argument_value_type;
            typedef const argument_value_type& argument_type;
            static PyObject* convert(argument_type val)
            {
                PyObject* retval =
                    PyTuple_New(boost::tuples::length<Ttuple_>::value);
                build_pytuple_from_tuple(retval, val);
                return retval;
            }
        };

        template<typename Tfirst_, typename Tsecond_>
        struct tuple_to_tuple_converter<std::pair<Tfirst_, Tsecond_> >
        {
            typedef std::pair<Tfirst_, Tsecond_> argument_value_type;
            typedef const argument_value_type& argument_type;

            static PyObject* convert(argument_type val)
            {
                return boost::python::incref(
                        boost::python::make_tuple(
                                val.first, val.second).ptr());
            }
        };

        template<typename Ttuple_>
        struct tuple_to_tuple_converter<boost::shared_ptr<Ttuple_> >
        {
            typedef Ttuple_ argument_value_type;
            typedef boost::shared_ptr<argument_value_type> argument_type;
            static PyObject* convert(argument_type val)
            {
                return tuple_to_tuple_converter<argument_value_type>::convert(
                        *val);
            }
        };

        template<typename T_>
        struct to_python_converter_fun
            : public std::unary_function<T_, boost::python::object>
        {
            typedef T_ argument_type;
            typedef boost::python::object result_type;

            result_type operator()(const argument_type& src)
            {
                return boost::python::object(src);
            }
        };

        template<typename Talloc_>
        class default_initializer_fun
            : public std::unary_function<typename Talloc_::reference, void>
        {
        public:
            typedef typename Talloc_::reference argument_type;
            typedef void result_type;

        public:
            default_initializer_fun(Talloc_& alloc): alloc_(alloc) {}

            void operator()(argument_type ptr)
            {
                new(alloc_.address(ptr)) typename Talloc_::value_type();
            }

        private:
            Talloc_& alloc_;
        };

        template<typename Talloc_>
        class destructor_fun
            : public std::unary_function<typename Talloc_::reference, void>
        {
        public:
            typedef typename Talloc_::reference argument_type;
            typedef void result_type;

        public:
            destructor_fun(Talloc_& alloc): alloc_(alloc) {}

            void operator()(argument_type ptr)
            {
                typedef typename Talloc_::value_type value_type;
                ptr.~value_type();
            }

        private:
            Talloc_& alloc_;
        };

        namespace for_compile_time_error {
            template<typename T_>
            class numpy_does_not_support_the_type_;

            template<typename T_>
            class numpy_type_does_not_yield_from_;
        };

        template<typename T_>
        struct get_numpy_typecode {
            static const std::size_t _ = sizeof(
                    for_compile_time_error::
                    numpy_does_not_support_the_type_<T_>);
        };

#       define DEFINE_NUMPY_TYPECODE_ASSOC(__type__, __value__)  \
            template<> struct get_numpy_typecode<__type__> \
            { \
                BOOST_STATIC_CONSTANT(enum NPY_TYPES, value = __value__); \
            }
        DEFINE_NUMPY_TYPECODE_ASSOC(bool,            NPY_BOOL);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_byte,        NPY_BYTE);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_ubyte,       NPY_UBYTE);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_short,       NPY_SHORT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_ushort,      NPY_USHORT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_int,         NPY_INT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_uint,        NPY_UINT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_long,        NPY_LONG);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_ulong,       NPY_ULONG);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_longlong,    NPY_LONGLONG);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_ulonglong,   NPY_ULONGLONG);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_float,       NPY_FLOAT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_double,      NPY_DOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_longdouble,  NPY_LONGDOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_cfloat,      NPY_CFLOAT);
        DEFINE_NUMPY_TYPECODE_ASSOC(std::complex<npy_float>, NPY_CFLOAT);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_cdouble,     NPY_CDOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(std::complex<npy_double>, NPY_CDOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(npy_clongdouble, NPY_CLONGDOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(
            std::complex<npy_longdouble>, NPY_CLONGDOUBLE);
        DEFINE_NUMPY_TYPECODE_ASSOC(boost::python::object, NPY_OBJECT);
        DEFINE_NUMPY_TYPECODE_ASSOC(std::string,           NPY_STRING);
#       define TMP std::basic_string<wchar_t, std::char_traits<wchar_t> >
        DEFINE_NUMPY_TYPECODE_ASSOC(TMP, NPY_UNICODE);
#       undef TMP
        DEFINE_NUMPY_TYPECODE_ASSOC(void,                  NPY_VOID);
        DEFINE_NUMPY_TYPECODE_ASSOC(char,                  NPY_CHAR);
#       undef DEFINE_NUMPY_TYPECODE_ASSOC

        template<typename Tarray_>
        struct to_ndarray_converter
        {
            typedef for_compile_time_error::
                    numpy_type_does_not_yield_from_<Tarray_> _;
        };

        template<typename T_, std::size_t N_>
        struct to_ndarray_converter<
                boost::multi_array<T_, N_, pyarray_backed_allocator<T_> > >
        {
            typedef boost::multi_array<T_, N_, pyarray_backed_allocator<T_> > source_type;
            static PyObject* convert(const source_type& val)
            {
                const npy_intp* dims;
                npy_intp _dims[N_];

                if (sizeof(npy_intp) == sizeof(typename source_type::size_type))
                {
                    dims = reinterpret_cast<const npy_intp*>(val.shape());
                }
                else
                {
                    for (std::size_t i = 0; i < N_; ++i)
                        _dims[i] = val.shape()[i];
                    dims = _dims;
                }
                return PyArray_New(&PyArray_Type, N_,
                        const_cast<npy_intp*>(dims),
                        get_numpy_typecode<T_>::value, NULL,
                        const_cast<source_type&>(val).origin(), 0,
                        NPY_CARRAY | NPY_OWNDATA, NULL);
            }
        };

        template<typename T_>
        struct to_ndarray_converter<std::vector<T_, pyarray_backed_allocator<T_ > > >
        {
            typedef std::vector<T_, pyarray_backed_allocator<T_> > source_type;

            static PyObject* convert(const source_type& val)
            {
                const npy_intp dims[1] = { val.size() };

                return PyArray_New(&PyArray_Type, 1,
                        const_cast<npy_intp*>(dims),
                        get_numpy_typecode<T_>::value, NULL,
                        const_cast<source_type&>(val).data(), 0,
                        NPY_CARRAY | NPY_OWNDATA, NULL);
            }
        };

        template<typename T_>
        struct to_ndarray_converter<std::vector<T_> >
        {
            typedef std::vector<T_> source_type;

            static PyObject* convert(const source_type& val)
            {
                typedef pyarray_backed_allocator<boost::python::object>
                        pyobject_array_allocator_type;

                BOOST_STATIC_ASSERT(
                        sizeof(boost::python::object) == sizeof(PyObject*)); 

                const npy_intp dims[1] = { val.size() };
                pyobject_array_allocator_type alloc(false);
                boost::python::object* converted_data =
                        alloc.allocate(val.size());

                boost::python::object* di = converted_data; 
                try
                {
                    for (typename source_type::const_iterator i(val.begin()),
                                                              e(val.end());
                            i != e; ++i)
                    {
                        new(di) boost::python::object(*i);
                        ++di; // this must be incremented after the pointed
                              // object pointed is successfully initialized
                    }
                }
                catch (const std::exception&)
                {
                    std::for_each(converted_data, di, destructor_fun<
                            pyobject_array_allocator_type>(alloc));
                    return NULL; 
                }

                return PyArray_New(&PyArray_Type, 1,
                        const_cast<npy_intp*>(dims),
                        get_numpy_typecode<boost::python::object>::value, NULL,
                        converted_data, 0,
                        NPY_CARRAY | NPY_OWNDATA, NULL);
            }
        };
    }

    void register_stop_iteration_exc_translator()
    {
        boost::python::register_exception_translator<StopIteration>(
                &detail::_handle_stop_iteration_exc);
    }

    template<typename Ttuple_>
    void register_tuple_converter()
    {
        boost::python::to_python_converter<
            Ttuple_, detail::tuple_to_tuple_converter<Ttuple_> >();
        boost::python::to_python_converter<
            boost::shared_ptr<Ttuple_>,
            detail::tuple_to_tuple_converter<boost::shared_ptr<Ttuple_> > >();
    }

    template<typename Tarray_>
    void register_multi_array_converter()
    {
        boost::python::to_python_converter<
            Tarray_, detail::to_ndarray_converter<Tarray_> >();
    }

    template<typename Tnative_, typename Tconverter_>
    void to_native_converter()
    {
        boost::python::converter::registry::push_back(
                &Tconverter_::convertible,
                reinterpret_cast<
                        boost::python::converter::constructor_function>(
                            &Tconverter_::construct),
                boost::python::type_id<Tnative_>());
    }
/*
    template<typename Tgen_>
    Enumerator<typename boost::remove_pointer<
            typename Tgen_::result_type>::type&>*
    make_enumerator(const Tgen_& gen)
    {
        return new GeneratorEnumeratorAdapter<Tgen_>(gen);
    }
*/
} // namespace util


class Sphere
{
public:
    typedef ::sphere<double> impl_type;
    typedef impl_type::value_type value_type;

public:
    Sphere(): impl_(impl_type::position_type(0, 0, 0), 0) {}

    Sphere(const impl_type::position_type& p, value_type r): impl_(p, r) {}

    PyObject* __repr__() const
    {
        return util::pystr_from_repr(&impl_);
    }

    value_type getX() const
    {
        return impl_.position.x();
    }

    void setX(value_type val)
    {
        impl_.position.x() = val;
    }

    value_type getY() const
    {
        return impl_.position.y();
    }

    void setY(value_type val)
    {
        impl_.position.y() = val;
    }

    value_type getZ() const
    {
        return impl_.position.z();
    }

    void setZ(value_type val)
    {
        impl_.position.z() = val;
    }

    value_type getRadius() const
    {
        return impl_.radius;
    }

    void setRadius(value_type val)
    {
        impl_.radius = val;
    }

    operator impl_type&()
    {
        return impl_;
    }

    operator const impl_type&() const
    {
        return impl_;
    }

    struct PositionToNDArrayConverter
    {
        typedef impl_type::position_type native_type;

        static PyObject* convert(const native_type& p)
        {
            static const npy_intp dims[1] = { native_type::size() };
            return PyArray_New(&PyArray_Type, 1, const_cast<npy_intp*>(dims),
                   util::detail::get_numpy_typecode<value_type>::value, NULL,
                    const_cast<void*>(static_cast<const void*>(&p)),
                    0, NPY_CARRAY, NULL);
        }
    };

    struct NDArrayToPositionConverter
    {
        typedef impl_type::position_type native_type;

        static void* convertible(PyObject* ptr)
        {
            if (!PyArray_Check(ptr))
            {
                return NULL;
            }

            PyObject* retval(PyArray_CastToType(
                    reinterpret_cast<PyArrayObject*>(ptr),
                    PyArray_DescrFromType(
                        util::detail::get_numpy_typecode<
                            native_type::value_type>::value), 0));
            if (!retval)
            {
                return NULL;
            }

            if (PyObject_Size(retval) != 3)
            {
                boost::python::decref(retval);
                return NULL;
            }

            return retval;
        }

        static void construct(PyObject* ptr,
                boost::python::converter::rvalue_from_python_storage<native_type>* data)
        {
            PyArrayObject* array_obj = static_cast<PyArrayObject*>(
                    data->stage1.convertible);
            data->stage1.convertible = new(data->storage.bytes) native_type(
                    reinterpret_cast<double*>(PyArray_DATA(array_obj)));
            boost::python::decref(reinterpret_cast<PyObject*>(array_obj));
        }
    };

    struct SequenceToPositionConverter
    {
        typedef impl_type::position_type native_type;

        static void* convertible(PyObject* ptr)
        {
            if (!PySequence_Check(ptr))
            {
                return NULL;
            }

            if (PySequence_Size(ptr) != 3)
            {
                return NULL;
            }

            return ptr;
        }

        static void construct(PyObject* ptr,
                boost::python::converter::rvalue_from_python_storage<native_type>* data)
        {
            data->stage1.convertible = new(data->storage.bytes) native_type(
                PyFloat_AsDouble( PySequence_GetItem(ptr, 0) ),
                PyFloat_AsDouble( PySequence_GetItem(ptr, 1) ),
                PyFloat_AsDouble( PySequence_GetItem(ptr, 2) ) );

/*
            boost::array<PyObject*, 3> items = {
                PySequence_GetItem(ptr, 0),
                PySequence_GetItem(ptr, 1),
                PySequence_GetItem(ptr, 2)
            };

            data->stage1.convertible = new(data->storage.bytes) native_type(
                PyFloat_AsDouble(items[0]),
                PyFloat_AsDouble(items[1]),
                PyFloat_AsDouble(items[2]));
*/
        }
    };

    inline static void __register_class()
    {
        using namespace boost::python;

        to_python_converter<impl_type::position_type,
                PositionToNDArrayConverter>();
        util::to_native_converter<impl_type::position_type,
                NDArrayToPositionConverter>();
        util::to_native_converter<impl_type::position_type,
                SequenceToPositionConverter>();

        class_<Sphere>("Sphere")
            .def(init<const ::position<double>&, double>())
            .def("__repr__", &Sphere::__repr__)
            .add_property("x", &Sphere::getX, &Sphere::setX)
            .add_property("y", &Sphere::getY, &Sphere::setY)
            .add_property("z", &Sphere::getZ, &Sphere::setZ)
            .add_property("radius", &Sphere::getRadius, &Sphere::setRadius);
    }
private:
    impl_type impl_;
};

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

/*
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
*/

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
                da_.push_back(pos_.distance((*i).second.position));
            }

            inline void operator()(impl_type::iterator i,
                    const position_type& d)
            {
                sa_.push_back(i);
                da_.push_back(pos_.distance((*i).second.position + d));
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

        value_type getX() const
        {
            return (*impl_).second.position.x();
        }

        void setX(value_type val)
        {
            (*impl_).second.position.x() = val;
        }

        value_type getY() const
        {
            return (*impl_).second.position.y();
        }

        void setY(value_type val)
        {
            (*impl_).second.position.y() = val;
        }

        value_type getZ() const
        {
            return (*impl_).second.position.z();
        }

        void setZ(value_type val)
        {
            (*impl_).second.position.z() = val;
        }

        value_type getRadius() const
        {
            return (*impl_).second.radius;
        }

        void setRadius(value_type val)
        {
            (*impl_).second.radius = val;
        }

        boost::python::object getID() const
        {
            return (*impl_).first;
        }
   
        inline static void __register_class()
        {
            using namespace boost::python;

            class_<SphereRef>("SphereRef", no_init)
                .def("__repr__", &SphereRef::__repr__)
                .add_property("x", &SphereRef::getX, &SphereRef::setX)
                .add_property("y", &SphereRef::getY, &SphereRef::setY)
                .add_property("z", &SphereRef::getZ, &SphereRef::setZ)
                .add_property("radius",
                        &SphereRef::getRadius,
                        &SphereRef::setRadius)
                .add_property("id", &SphereRef::getID);
        }
 
    private:
        impl_type impl_;
    };

    struct IterToSphereRefConverter
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

/*
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
*/

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
                IterToSphereRefConverter>();

//        util::register_tuple_converter<Generators::result_type>();

        util::register_multi_array_converter<
            boost::tuples::element<0, Builders::result_type>::type>();
        util::register_multi_array_converter<
            boost::tuples::element<1, Builders::result_type>::type>();
        // the following conversion is the same as the previous
        // util::register_multi_array_converter<boost::tuples::element<2, Builders::result_type>::type>();

        util::register_tuple_converter<Builders::result_type>();
/*
        util::register_enumerator<
                const Generators::result_type&,
                return_value_policy<copy_const_reference> >(
                "ObjectContainer_NeighborIterator");
*/
        class_<ObjectContainer>("ObjectContainer")
            .def(init<length_type, matrix_size_type>())
            .add_property("cell_size", &ObjectContainer::cell_size)
            .add_property("world_size", &ObjectContainer::world_size)
            .add_property("matrix_size", &ObjectContainer::matrix_size)
/*
            .def("iterneighbors", &ObjectContainer::iterneighbors,
                    return_value_policy<manage_new_object>())
            .def("iterneighbors_cyclic", &ObjectContainer::iterneighbors_cyclic,
                    return_value_policy<manage_new_object>())
*/
            .def("neighbors_array", &ObjectContainer::neighbors_array)
            .def("neighbors_array_cyclic", &ObjectContainer::neighbors_array_cyclic)
            .def("all_neighbors_array", &ObjectContainer::all_neighbors_array)
            .def("all_neighbors_array_cyclic", &ObjectContainer::all_neighbors_array_cyclic)
            .def("__len__", &ObjectContainer::__len__)
            .def("__setitem__", &ObjectContainer::__setitem__)
            .def("__getitem__", &ObjectContainer::__getitem__,
                        return_value_policy<manage_new_object>());
    }

private:
    impl_type impl_;
};

} // namespace peer


BOOST_PYTHON_MODULE(_object_matrix)
{
    using namespace boost::python;

    import_array();
    peer::util::register_stop_iteration_exc_translator();
    peer::Sphere::__register_class();
    peer::ObjectContainer::__register_class();
}
