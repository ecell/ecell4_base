#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "pointer_as_ref"

#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>

#include "utils/pointer_as_ref.hpp"

BOOST_AUTO_TEST_CASE(basic)
{
    int a = 0;
    pointer_as_ref<int, int*> b(&a);
    static_cast<int&>(b) = 1;
    BOOST_CHECK_EQUAL(a, 1);
}

BOOST_AUTO_TEST_CASE(reference_holder)
{
    int a = 0;
    int *b = 0;
    pointer_as_ref<int, int*&> c(b);
    b = &a;
    static_cast<int&>(c) = 1;
    BOOST_CHECK_EQUAL(a, 1);
}
