#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "range_support_test"

#include <boost/test/included/unit_test.hpp>
#include <list>
#include <vector>
#include "utils/range.hpp"
#include "utils/range_support.hpp"
    
struct fn_t {
    void operator()(int v)
    {
        value = v;
    }

    fn_t(): value(-1) {}

    int value;
};

BOOST_AUTO_TEST_CASE(test_check_range_iterator_category_mf)
{
    BOOST_CHECK((check_range_iterator_category<int[4], boost::random_access_traversal_tag>::value));
    BOOST_CHECK((check_range_iterator_category<std::vector<int>, boost::random_access_traversal_tag>::value));
    BOOST_CHECK((!check_range_iterator_category<std::list<int>, boost::random_access_traversal_tag>::value));
}

BOOST_AUTO_TEST_CASE(test_call_with_size_if_randomly_accessible)
{
    {
        int test[4];
        fn_t fn;
        call_with_size_if_randomly_accessible(fn, test);
        BOOST_CHECK_EQUAL(fn.value, 4);
    }

    {
        fn_t fn;
        call_with_size_if_randomly_accessible(fn, std::list<int>());
        BOOST_CHECK_EQUAL(fn.value, -1);
    }
}

