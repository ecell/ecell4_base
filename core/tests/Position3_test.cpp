#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>

#include "../Position3.hpp"
#include "../linear_algebra.hpp"

using namespace ecell4;


class Position3Test
    : public CppUnit::TestFixture
{
public:

    CPPUNIT_TEST_SUITE(Position3Test);
    CPPUNIT_TEST(test_constructor);
    CPPUNIT_TEST(test_four_arithmetic_operations);
    CPPUNIT_TEST_SUITE_END();

public:

    void setUp()
    {
        target = new Position3();
    }

    void tearDown()
    {
        delete target;
    }

    void test_constructor();
    void test_four_arithmetic_operations();
    void test_add();
    void test_substract();

private:

    Position3 *target;
};

CPPUNIT_TEST_SUITE_REGISTRATION(Position3Test);

void Position3Test::test_constructor()
{
    CPPUNIT_ASSERT_EQUAL(*target, Position3(0, 0, 0));
}

void Position3Test::test_four_arithmetic_operations()
{
    Position3 pos1(1, 2, 3);
    CPPUNIT_ASSERT_EQUAL(pos1 * 2, Position3(2, 4, 6));
}

void Position3Test::test_add()
{
    Position3 pos1(1,2,3);
    Position3 pos2(2,4,6);
    CPPUNIT_ASSERT_EQUAL(pos1 + pos2, Position3(3, 6, 9));
}

void Position3Test::test_substract()
{
    Position3 pos2(2,4,6);
    Position3 pos1(1,2,3);
    CPPUNIT_ASSERT_EQUAL(pos2 - pos1, Position3(1, 2, 3));
}


