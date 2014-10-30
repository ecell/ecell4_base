#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "py_range_converters_test"

#include <boost/test/included/unit_test.hpp>
#include "peer/wrappers/range/pyiterable_range.hpp"

BOOST_AUTO_TEST_CASE(pyiterable_range_test)
{
    Py_InitializeEx(0);
    {
        boost::python::list pylist;
        pylist.append(1);
        pylist.append(2);
        pylist.append(3);
        peer::wrappers::pyiterable_range<int> wrapper(pylist);

        BOOST_CHECK_EQUAL(size(wrapper), boost::python::len(pylist));
    }
    Py_Finalize();
}

