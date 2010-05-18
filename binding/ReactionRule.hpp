#ifndef PEER_REACTION_RULE_HPP
#define PEER_REACTION_RULE_HPP

#include <iterator>
#include <algorithm>
#include <boost/python.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include "peer/converters/tuple.hpp"
#include "peer/converters/sequence.hpp"
#include "../ReactionRule.hpp"

namespace binding {

template<typename Trr_>
struct seq_to_reactants_converter
{
    typedef typename Trr_::Reactants native_type;
    
    static void* convertible(PyObject* ptr)
    {
        if (!PySequence_Check(ptr))
        {
            return NULL;
        }

        Py_ssize_t nitems = PySequence_Size(ptr);
        if (nitems < 1 || nitems > 2)
        {
            return NULL;
        }
        
        return ptr;
    }
    
    static void construct(PyObject* ptr,
                          boost::python::converter::rvalue_from_python_storage<native_type>* data)
    {
        using namespace boost::python;
        Py_ssize_t nitems = PySequence_Size(ptr);
        switch (nitems)
        {
        case 1:
            data->stage1.convertible = new(data->storage.bytes) native_type(
                extract<typename native_type::value_type>(
                    object(handle<>(PySequence_GetItem(ptr, 0)))));
            break;
        case 2:
            data->stage1.convertible = new(data->storage.bytes) native_type(
                extract<typename native_type::value_type>(
                    object(handle<>(PySequence_GetItem(ptr, 0)))),
                extract<typename native_type::value_type>(
                    object(handle<>(PySequence_GetItem(ptr, 1)))));
            break;
        default:
            // never get here 
            break;
        }
    }
};

template<typename Timpl_>
struct ReactionRuleExtras
{
    typedef Timpl_ impl_type;

    static std::string __str__(impl_type* impl)
    {
        return boost::lexical_cast<std::string>(*impl);
    }

    static bool __eq__(PyObject* self, PyObject* rhs)
    {
        using namespace boost::python;
        impl_type *impl = extract<impl_type*>(object(borrowed(self)));
        if (!PyObject_TypeCheck(rhs, self->ob_type))
        {
            return false;
        }
        return *impl == *extract<impl_type*>(object(borrowed(rhs)));
    }

    static typename impl_type::species_type_id_iterator get_products_begin(impl_type *impl)
    {
        return impl->get_products().begin();
    }

    static typename impl_type::species_type_id_iterator get_products_end(impl_type *impl)
    {
        return impl->get_products().end();
    }

    static std::string const& __getitem__(impl_type* impl, std::string const& key)
    {
        return (*impl)[key];
    }

    static void __setitem__(impl_type* impl, std::string const& key, std::string const& val)
    {
        (*impl)[key] = val;
    }
};

template<typename Trr_>
inline void register_reaction_rule_class(const char* name)
{
    using namespace boost::python;
    typedef Trr_ impl_type;
    typedef ReactionRuleExtras<impl_type> extras_type;

    typedef std::vector<typename impl_type::species_type_id_type> species_type_id_vector;
    peer::util::to_native_converter<
            typename impl_type::Reactants,
            seq_to_reactants_converter<impl_type> >();
    peer::converters::register_range_to_tuple_converter<
            typename impl_type::Reactants>();
    peer::converters::register_iterable_to_range_converter<
            species_type_id_vector>();

    class_<impl_type, impl_type*>(name,
            init<typename impl_type::Reactants const&, species_type_id_vector>())
        .add_property("reactants",  
            make_function(&impl_type::get_reactants,
                return_value_policy<copy_const_reference>()))
        .add_property("products",
            range<return_value_policy<return_by_value>, impl_type*>(
                  &extras_type::get_products_begin,
                  &extras_type::get_products_end))
        .def("__getitem__", &extras_type::__getitem__,
                return_value_policy<copy_const_reference>())
        .def("__setitem__", &extras_type::__setitem__)
        .def("__str__", &extras_type::__str__)
        .def("__eq__", &extras_type::__eq__)
        ;
}

} // namespace binding

#endif /* PEER_REACTION_RULE_HPP */
