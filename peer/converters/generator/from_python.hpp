#ifndef PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP
#define PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP

namespace peer { namespace converters {

template<typename Twrapper_, typename Tgen_ = typename Twrapper_::generator_type>
struct converter_pair
{
    typedef Twrapper_ self_type;
    typedef Tgen_ generator_type;

    struct to_python_converter
    {
        static PyObject* convert(Tgen_ const& impl)
        {
            return self_type::create(impl);
        }

        static PyTypeObject* get_pytype()
        {
            return &self_type::__class__;
        }
    };

    struct to_native_converter
    {
        static void* convertible(PyObject* pyo)
        {
            if (!PyObject_TypeCheck(pyo, &self_type::__class__))
            {
                return 0;
            }
            return pyo;
        }

        static void construct(PyObject* pyo, 
                              boost::python::converter::rvalue_from_python_stage1_data* data)
        {
            void* storage(reinterpret_cast<
                boost::python::converter::rvalue_from_python_storage<Tgen_>* >(
                    data)->storage.bytes);
            new (storage) Tgen_(reinterpret_cast<self_type*>(pyo)->impl_);
            data->convertible = storage;
        }
    };

    static void __register_converter()
    {
        util::to_native_converter<Tgen_, to_native_converter>();
        boost::python::to_python_converter<Tgen_, to_python_converter>();
    }
};

template<typename Twrapper_, typename Tgen_, typename Tholder_>
struct converter_pair<Twrapper_, ptr_generator<Tgen_, Tholder_> >
{
    typedef Twrapper_ self_type;
    typedef ptr_generator<Tgen_, Tholder_> generator_type;

    struct to_python_converter
    {
        static PyObject* convert(Tgen_* impl)
        {
            if (impl)
            {
                Tholder_ ptr(impl);
                return self_type::create(generator_type(ptr));
            }
            return boost::python::incref(Py_None);
        }

        static PyTypeObject* get_pytype()
        {
            return &self_type::__class__;
        }
    };

    struct to_native_converter
    {
        static void* convert(PyObject* pyo)
        {
            return reinterpret_cast<self_type*>(pyo)->impl_.ptr().get();
        }

        static PyTypeObject const* expected_pytype()
        {
            return &self_type::__class__;
        }
    };

    static void __register_converter()
    {
        util::to_native_lvalue_converter<Tgen_*, to_native_converter>();
        boost::python::to_python_converter<Tgen_*, to_python_converter>();
    }
};

} // namespace peer::converters

#endif /* PEER_CONVERTERS_GENERATOR_TO_PYTHON_HPP */
