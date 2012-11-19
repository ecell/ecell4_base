#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include "../BDWorld.hpp"

using namespace ecell4;
using namespace ecell4::bd;


class BDWorldTest
    : public CppUnit::TestFixture
{
public:

    CPPUNIT_TEST_SUITE(BDWorldTest);
    CPPUNIT_TEST(test_edge_lengths);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp()
    {
        Real const L(1e-6);
        Position3 const edge_lengths(L, L, L);
        target = new BDWorld(edge_lengths);
    }

    void tearDown()
    {
        delete target;
    }

    void test_edge_lengths();

private:

    BDWorld *target;
};

CPPUNIT_TEST_SUITE_REGISTRATION(BDWorldTest);

void BDWorldTest::test_edge_lengths()
{
    Position3 const& edge_lengths(target->edge_lengths());
    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        CPPUNIT_ASSERT(edge_lengths[dim] > 0);
    }
}
