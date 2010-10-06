#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "binding_common.hpp"
#include "StructureType.hpp"
#include "peer/utils.hpp"

namespace binding {

static boost::python::object structure_type_class;

template<typename Tsid_, typename Tst_>
struct structure_type_to_structure_id_converter
{
    typedef Tst_ structure_type_type;
    typedef Tsid_ native_type;

    static void* convertible(PyObject* pyo)
    {
        if (!PyObject_TypeCheck(pyo, reinterpret_cast<PyTypeObject*>(
                structure_type_class.ptr())))
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
        new (storage) native_type(extract<structure_type_type const*>(object(borrowed(pyo)))()->id());
        data->convertible = storage;
    }
};

void register_structure_type_class()
{
    structure_type_class = register_structure_type_class<StructureType>("StructureType");
    peer::util::to_native_converter<StructureID, structure_type_to_structure_id_converter<StructureID, StructureType> >();
}

} // namespace binding
