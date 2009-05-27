#define BOOST_TEST_MODULE "vmap"
#include <boost/test/included/unit_test.hpp>

#include "vmap.hpp"

BOOST_AUTO_TEST_CASE(insert)
{
    typedef vmap<int, int> map_type;
    map_type m;

    BOOST_CHECK_EQUAL(0, m.size());
    m.insert(std::make_pair(1, 0));
    BOOST_CHECK_EQUAL(1, m.size());
    BOOST_CHECK_EQUAL(0, *(m.values().end() - 1));
    m.insert(std::make_pair(2, 1));
    BOOST_CHECK_EQUAL(2, m.size());
    BOOST_CHECK_EQUAL(1, *(m.values().end() - 1));
    m.insert(std::make_pair(3, 2));
    BOOST_CHECK_EQUAL(3, m.size());
    BOOST_CHECK_EQUAL(2, *(m.values().end() - 1));
    m.insert(std::make_pair(4, 3));
    BOOST_CHECK_EQUAL(4, m.size());
    BOOST_CHECK_EQUAL(3, *(m.values().end() - 1));
}

BOOST_AUTO_TEST_CASE(subscription)
{
    typedef vmap<int, int> map_type;
    map_type m;

    BOOST_CHECK_EQUAL(0, m.size());
    m.insert(std::make_pair(1, 0));
    BOOST_CHECK_EQUAL(1, m.size());
    BOOST_CHECK_EQUAL(*(m.values().begin() + 0), m[1]);
    m.insert(std::make_pair(2, 1));
    BOOST_CHECK_EQUAL(2, m.size());
    BOOST_CHECK_EQUAL(*(m.values().begin() + 1), m[2]);
    m.insert(std::make_pair(3, 2));
    BOOST_CHECK_EQUAL(3, m.size());
    BOOST_CHECK_EQUAL(*(m.values().begin() + 2), m[3]);
    m.insert(std::make_pair(4, 3));
    BOOST_CHECK_EQUAL(4, m.size());
    BOOST_CHECK_EQUAL(*(m.values().begin() + 3), m[4]);
}

BOOST_AUTO_TEST_CASE(find)
{
    typedef vmap<int, int> map_type;
    map_type m;

    BOOST_CHECK_EQUAL(0, m.size());
    BOOST_CHECK(m.end() == m.find(0));
    BOOST_CHECK(m.end() == m.find(1));
    BOOST_CHECK(m.end() == m.find(2));
    BOOST_CHECK(m.end() == m.find(3));
    BOOST_CHECK(m.end() == m.find(4));
    m.insert(std::make_pair(1, 0));
    BOOST_CHECK_EQUAL(1, m.size());
    BOOST_CHECK(m.end() == m.find(0));
    BOOST_CHECK(m.end() != m.find(1));
    BOOST_CHECK_EQUAL(0, (*m.find(1)).second);
    m.insert(std::make_pair(2, 1));
    BOOST_CHECK_EQUAL(2, m.size());
    BOOST_CHECK(m.end() == m.find(0));
    BOOST_CHECK(m.end() != m.find(1));
    BOOST_CHECK(m.end() != m.find(2));
    BOOST_CHECK_EQUAL(1, (*m.find(2)).second);
    m.insert(std::make_pair(3, 2));
    BOOST_CHECK_EQUAL(3, m.size());
    BOOST_CHECK(m.end() == m.find(0));
    BOOST_CHECK(m.end() != m.find(1));
    BOOST_CHECK(m.end() != m.find(2));
    BOOST_CHECK(m.end() != m.find(3));
    BOOST_CHECK_EQUAL(2, (*m.find(3)).second);
    m.insert(std::make_pair(4, 3));
    BOOST_CHECK_EQUAL(4, m.size());
    BOOST_CHECK(m.end() == m.find(0));
    BOOST_CHECK(m.end() != m.find(1));
    BOOST_CHECK(m.end() != m.find(2));
    BOOST_CHECK(m.end() != m.find(3));
    BOOST_CHECK(m.end() != m.find(4));
    BOOST_CHECK_EQUAL(3, (*m.find(4)).second);
}

BOOST_AUTO_TEST_CASE(erase)
{
    typedef vmap<int, int> map_type;
    map_type m;

    m.insert(std::make_pair(1, 0));
    m.insert(std::make_pair(2, 1));
    m.insert(std::make_pair(3, 2));
    m.insert(std::make_pair(4, 3));
    BOOST_CHECK_EQUAL(4, m.size());
    m.erase(m.find(1));
    BOOST_CHECK_EQUAL(3, m.size());
    BOOST_CHECK(m.end() == m.find(0));
    BOOST_CHECK(m.end() == m.find(1));
    BOOST_CHECK(m.end() != m.find(2));
    BOOST_CHECK(m.end() != m.find(3));
    BOOST_CHECK(m.end() != m.find(4));
}

