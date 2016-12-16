#include "PyapiBind.hpp"

namespace ecell4
{

PythonAPIBind::PythonAPIBind()
{
    ;
}

PythonAPIBind::PythonAPIBind(PyRefFuncType inc_ref, PyRefFuncType dec_ref):
    inc_ref_(inc_ref), dec_ref_(dec_ref)
{
    ;
}


} // namespace ecell4
