#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/python.hpp>
#include "peer/wrappers/range/stl_container_wrapper.hpp"
#include "binding_common.hpp"
#include "NetworkRulesWrapper.hpp"
#include "ReactionRuleInfo.hpp"
#include "SpeciesInfo.hpp"

namespace binding {

static boost::python::object species_info_class;

template<typename Timpl_>
struct reaction_rule_vector_converter
{
    typedef Timpl_ native_type;

    struct instance_holder
    {
        instance_holder(native_type const& instance): instance_(instance) {}

        native_type const& operator*() const
        {
            return instance_;
        }

        native_type const* operator->() const
        {
            return &(**this);
        }

        native_type& operator*()
        {
            PyErr_SetString(PyExc_RuntimeError, "object is immutable");
            boost::python::throw_error_already_set();
            return *static_cast<native_type*>(0);
        }

        native_type* operator->()
        {
            return &(**this);
        }

        native_type const& instance_;
    };

    typedef peer::wrappers::stl_container_wrapper<native_type, instance_holder> wrapper_type;

    struct to_python_converter
    {
        static PyObject* convert(native_type const& v)
        {
            return reinterpret_cast<PyObject*>(wrapper_type::create(instance_holder(v)));
        }
    };

    struct to_native_converter
    {
        static void* convert(PyObject* ptr)
        {
            return const_cast<native_type*>(&*reinterpret_cast<wrapper_type const*>(ptr)->ptr());
        }

        static PyTypeObject const* expected_pytype()
        {
            return &wrapper_type::__class__;
        }
    };

    static void __register()
    {
        wrapper_type::__class_init__("ReactionRuleInfoVector", boost::python::scope().ptr());
        boost::python::to_python_converter<native_type, to_python_converter>();
        peer::util::to_native_lvalue_converter<native_type, to_native_converter>();
    }
};

template<typename Tsid_, typename Tsinfo_>
struct species_info_to_species_id_converter
{
    typedef Tsinfo_ species_info_type;
    typedef Tsid_ native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyObject_TypeCheck(pyo, reinterpret_cast<PyTypeObject*>(
                species_info_class.ptr())))
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
            converter::rvalue_from_python_storage<native_type>* >(
                data)->storage.bytes);
        new (storage) native_type(static_cast<species_info_type*>(extract<species_info_type*>(object(borrowed(pyo))))->id());
        data->convertible = storage;
    }
};

static void register_reaction_rule_info_class()
{
    register_reaction_rule_info_class<ReactionRuleInfo>("ReactionRuleInfo");
}

static void register_species_info_class()
{
    species_info_class = register_species_info_class<SpeciesInfo>("SpeciesInfo");
    peer::util::to_native_converter<SpeciesID, species_info_to_species_id_converter<SpeciesID, SpeciesInfo> >();
}

void register_network_rules_wrapper_class()
{
    using namespace boost::python;

    register_network_rules_wrapper_class<NetworkRulesWrapper>("NetworkRulesWrapper");
    register_reaction_rule_info_class();
    reaction_rule_vector_converter<NetworkRulesWrapper::reaction_rules>::__register();
    register_species_info_class();
}

} // namespace binding
