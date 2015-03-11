#define BOOST_TEST_MODULE "Shape_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <boost/test/floating_point_comparison.hpp>

#include <ecell4/core/Sphere.hpp>
#include <ecell4/core/Rod.hpp>

using namespace ecell4;

struct Fixture
{
    const Real3 center;
    const Real radius;
    Sphere sphere;
    Fixture() :
        center(2.5e-6, 2.5e-6, 2.5e-6),
        radius(2.5e-7), sphere(center, radius)
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(suite, Fixture)

BOOST_AUTO_TEST_CASE(Shape_test_constructor)
{
}

BOOST_AUTO_TEST_CASE(Shape_test_is_inside)
{
    BOOST_CHECK(sphere.is_inside(center) <= 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.3e-6, 2.5e-6, 2.5e-6)) <= 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.5e-6, 2.3e-6, 2.5e-6)) <= 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.5e-6, 2.5e-6, 2.3e-6)) <= 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.2e-6, 2.5e-6, 2.5e-6)) > 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.5e-6, 2.2e-6, 2.5e-6)) > 0);
    BOOST_CHECK(sphere.is_inside(Real3(2.5e-6, 2.5e-6, 2.2e-6)) > 0);
}

BOOST_AUTO_TEST_SUITE_END()


struct RodFixture
{
    const Real3 center;
    const Real length;
    const Real radius;
    Rod rod;
    boost::shared_ptr<RandomNumberGenerator> rng;
    RodFixture() :
        center(5e-6, 5e-6, 5e-6), length(2.5e-6),
        radius(1.25e-6), rod(length, radius, center),
        rng(new GSLRandomNumberGenerator())
    {
    }
};

BOOST_FIXTURE_TEST_SUITE(rod_suite, RodFixture)

BOOST_AUTO_TEST_CASE(Rod_test_draw_position)
{
    for (int i(0); i < 1000; ++i)
        BOOST_ASSERT(rod.is_inside(rod.draw_position(rng)) <= 0);
}

BOOST_AUTO_TEST_CASE(RodSurface_test_draw_position)
{
    RodSurface surface(rod.surface());
    int over(0), just(0), under(0);
    for (int i(0); i < 1000; ++i)
    {
        const Real l(surface.is_inside(surface.draw_position(rng)));
        if (l > radius*1e-6)
            ++over;
        else if (l < -radius*1e-6)
            ++under;
        else
            ++just;
    }
    // std::cout << "over: " << over << ", just:" << just << ", under: " << under << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
