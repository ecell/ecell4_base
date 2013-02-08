#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include "../types.hpp"
#include "../ParticleSpace.hpp"

using namespace ecell4;


template<typename Timpl_>
class ParticleSpaceTest
    : public CppUnit::TestFixture
{
public:

    typedef Timpl_ implementation_type;

    CPPUNIT_TEST_SUITE(ParticleSpaceTest);
    CPPUNIT_TEST(test_edge_lengths);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp()
    {
        Real const L(1e-6);
        Position3 const edge_lengths(L, L, L);
        target = new implementation_type(edge_lengths);
    }

    void tearDown()
    {
        delete target;
    }

    void test_edge_lengths();

private:

    ParticleSpace *target;
};

CPPUNIT_TEST_SUITE_REGISTRATION(ParticleSpaceTest<ParticleSpaceVectorImpl>);

template<typename Timpl_>
void ParticleSpaceTest<Timpl_>::test_edge_lengths()
{
    Position3 const& edge_lengths(target->edge_lengths());
    for (Position3::size_type dim(0); dim < 3; ++dim)
    {
        CPPUNIT_ASSERT(edge_lengths[dim] > 0);
    }
}
