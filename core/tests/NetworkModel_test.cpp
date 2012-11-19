#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include "../types.hpp"
#include "../NetworkModel.hpp"

using namespace ecell4;


class NetworkModelTest
    : public CppUnit::TestFixture
{
public:

    CPPUNIT_TEST_SUITE(NetworkModelTest);
    CPPUNIT_TEST(test_add_species);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp()
    {
        target = new NetworkModel();
    }

    void tearDown()
    {
        delete target;
    }

    void test_add_species();

private:

    NetworkModel *target;
};

CPPUNIT_TEST_SUITE_REGISTRATION(NetworkModelTest);

void NetworkModelTest::test_add_species()
{
    ; // add a test case
}
