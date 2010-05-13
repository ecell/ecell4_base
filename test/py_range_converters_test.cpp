#define BOOST_TEST_MODULE "py_range_converters_test"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <boost/test/included/unit_test.hpp>
#include "peer/range_converters.hpp"

BOOST_AUTO_TEST_CASE(py_range_wrapper_test)
{
    Py_InitializeEx(0);
    {
        boost::python::list pylist;
        pylist.append(1);
        pylist.append(2);
        pylist.append(3);
        peer::util::py_range_wrapper<int> wrapper(pylist);

        BOOST_CHECK_EQUAL(boost::size(wrapper), boost::python::len(pylist));
    }
    Py_Finalize();
}

