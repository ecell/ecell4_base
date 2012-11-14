#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include "../types.hpp"
#include "../CompartmentSpace.hpp"

using namespace ecell4;


class CompartmentSpaceVectorImplTest
    : public CppUnit::TestFixture
{
public:

    CPPUNIT_TEST_SUITE(CompartmentSpaceVectorImplTest);
    CPPUNIT_TEST(test_volume);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp()
    {
        Real const volume(1e-18);
        target = new CompartmentSpaceVectorImpl(volume);
    }

    void tearDown()
    {
        delete target;
    }

    void test_volume();

private:

    CompartmentSpaceVectorImpl *target;
};

CPPUNIT_TEST_SUITE_REGISTRATION(CompartmentSpaceVectorImplTest);

void CompartmentSpaceVectorImplTest::test_volume()
{
    target->set_volume(2 * target->volume());
}
