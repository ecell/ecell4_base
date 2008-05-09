#include "peer/Sphere.hpp"
#include "peer/ObjectContainer.hpp"

BOOST_PYTHON_MODULE(object_matrix)
{
    using namespace boost::python;

    import_array();
#if OBJECTMATRIX_USE_ITERATOR
    peer::util::register_stop_iteration_exc_translator();
#endif
    peer::Sphere::__register_class();
    peer::ObjectContainer::__register_class();
}
