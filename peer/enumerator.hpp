#ifndef OBJECTMATRIX_PEER_ENUMERATOR_HPP
#define OBJECTMATRIX_PEER_ENUMERATOR_HPP

#include <stdexcept>

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
    namespace detail
    {
        inline void _handle_stop_iteration_exc(const StopIteration& exc)
        {
            PyErr_SetNone(PyExc_StopIteration);
        }
    } // namespace detail

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

    template<typename T_, typename Tpol_>
    void register_enumerator(const char* python_name)
    {
        boost::python::class_<EnumeratorWrapper<T_> >(python_name,
            boost::python::no_init)
            .def("next", &Enumerator<T_>::next, Tpol_())
            .def("__iter__", &pass_through);
    }

    void register_stop_iteration_exc_translator()
    {
        boost::python::register_exception_translator<StopIteration>(
                &detail::_handle_stop_iteration_exc);
    }

    template<typename Tgen_>
    Enumerator<typename boost::remove_pointer<
            typename Tgen_::result_type>::type&>*
    make_enumerator(const Tgen_& gen)
    {
        return new GeneratorEnumeratorAdapter<Tgen_>(gen);
    }
} // namespace util

} // namespace peer
#endif

#endif /* OBJECTMATRIX_PEER_ENUMERATOR_HPP */
