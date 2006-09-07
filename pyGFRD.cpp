#include <boost/python.hpp>

#include "PairGreensFunction.hpp"
#include "PlainPairGreensFunction.hpp"

using namespace boost::python;


BOOST_PYTHON_MODULE(gfrd)
{
  class_<PlainPairGreensFunction>( "PlainPairGreensFunction",
				   init<const Real, const Real, const Real>() )
    .def("getD", &PlainPairGreensFunction::getD)
    .def("getkf", &PlainPairGreensFunction::getkf)
    .def("getSigma", &PlainPairGreensFunction::getSigma)
    .def("drawR", &PlainPairGreensFunction::drawR )
    .def("drawTheta", &PlainPairGreensFunction::drawTheta )
    ;
}
