#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "binding_common.hpp"
#include "SpeciesType.hpp"
#include "peer/utils.hpp"

namespace binding {

static boost::python::object species_type_class;

template<typename Tsid_, typename Tst_>
struct species_type_to_species_id_converter
{
    typedef Tst_ species_type_type;
    typedef Tsid_ native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyObject_TypeCheck(pyo, reinterpret_cast<PyTypeObject*>(
                species_type_class.ptr())))
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
        new (storage) native_type(static_cast<species_type_type*>(extract<species_type_type*>(object(borrowed(pyo))))->id());
        data->convertible = storage;
    }
};

void register_species_type_class()
{
    species_type_class = register_species_type_class<SpeciesType>("SpeciesType");
    peer::util::to_native_converter<SpeciesID, species_type_to_species_id_converter<SpeciesID, SpeciesType> >();
}

} // namespace binding
