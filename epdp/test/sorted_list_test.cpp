#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "sorted_list"

#include <boost/test/included/unit_test.hpp>
#include <vector>
#include "sorted_list.hpp"

BOOST_AUTO_TEST_CASE(regression)
{
    sorted_list<std::vector<int> > a;
    a.push(299);
    a.push(300);
    a.push(301);
    BOOST_CHECK(a.erase(299));
    BOOST_CHECK(a.end() == a.find(299));
    BOOST_CHECK(a.end() != a.find(300));
    BOOST_CHECK(a.end() != a.find(301));
}
