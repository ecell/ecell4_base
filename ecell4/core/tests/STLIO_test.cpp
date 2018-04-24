#define BOOST_TEST_MODULE "STLIO_test"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#   include <boost/test/unit_test.hpp>
#else
#   define BOOST_TEST_NO_LIB
#   include <boost/test/included/unit_test.hpp>
#endif

#include <ecell4/core/STLFileIO.hpp>
#include <boost/random.hpp>
#include <utility>

using ecell4::Real;
using ecell4::Real3;
using ecell4::Triangle;

BOOST_AUTO_TEST_CASE(test_io_ascii)
{
    boost::random::mt19937 mt(123456789);
    boost::random::uniform_real_distribution<Real> uni(0.0, 100.0);

    std::vector<Triangle> triangles;
    for(std::size_t i=0; i<100; ++i)
    {
        const Real3 v0(uni(mt), uni(mt), uni(mt));
        const Real3 v1(uni(mt), uni(mt), uni(mt));
        const Real3 v2(uni(mt), uni(mt), uni(mt));
        triangles.push_back(Triangle(v0, v1, v2));
    }

    ecell4::write_stl_format(
            "STLIO_test_asc.stl", triangles, ecell4::STLFormat::Ascii);

    const std::vector<Triangle> after_io = ecell4::read_stl_format(
            "STLIO_test_asc.stl", ecell4::STLFormat::Ascii);

    BOOST_REQUIRE(after_io.size() == triangles.size());
    for(std::size_t i=0; i<triangles.size(); ++i)
    {
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(0)[0],
                                    after_io.at(i).vertex_at(0)[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(0)[1],
                                    after_io.at(i).vertex_at(0)[1], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(0)[1],
                                    after_io.at(i).vertex_at(0)[1], 1e-12);

        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(1)[0],
                                    after_io.at(i).vertex_at(1)[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(1)[1],
                                    after_io.at(i).vertex_at(1)[1], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(1)[2],
                                    after_io.at(i).vertex_at(1)[2], 1e-12);

        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(2)[0],
                                    after_io.at(i).vertex_at(2)[0], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(2)[1],
                                    after_io.at(i).vertex_at(2)[1], 1e-12);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(2)[2],
                                    after_io.at(i).vertex_at(2)[2], 1e-12);
    }
}

BOOST_AUTO_TEST_CASE(test_io_binary)
{
    boost::random::mt19937 mt(123456789);
    boost::random::uniform_real_distribution<Real> uni(0.0, 100.0);

    std::vector<Triangle> triangles;
    for(std::size_t i=0; i<100; ++i)
    {
        const Real3 v0(uni(mt), uni(mt), uni(mt));
        const Real3 v1(uni(mt), uni(mt), uni(mt));
        const Real3 v2(uni(mt), uni(mt), uni(mt));
        triangles.push_back(Triangle(v0, v1, v2));
    }

    ecell4::write_stl_format(
            "STLIO_test_bin.stl", triangles, ecell4::STLFormat::Binary);

    const std::vector<Triangle> after_io = ecell4::read_stl_format(
            "STLIO_test_bin.stl", ecell4::STLFormat::Binary);

    BOOST_REQUIRE(after_io.size() == triangles.size());
    for(std::size_t i=0; i<triangles.size(); ++i)
    {
        // because it uses float, the precision lost.
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(0)[0],
                                    after_io.at(i).vertex_at(0)[0], 1e-4);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(0)[1],
                                    after_io.at(i).vertex_at(0)[1], 1e-4);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(0)[1],
                                    after_io.at(i).vertex_at(0)[1], 1e-4);

        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(1)[0],
                                    after_io.at(i).vertex_at(1)[0], 1e-4);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(1)[1],
                                    after_io.at(i).vertex_at(1)[1], 1e-4);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(1)[2],
                                    after_io.at(i).vertex_at(1)[2], 1e-4);

        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(2)[0],
                                    after_io.at(i).vertex_at(2)[0], 1e-4);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(2)[1],
                                    after_io.at(i).vertex_at(2)[1], 1e-4);
        BOOST_CHECK_CLOSE_FRACTION(triangles.at(i).vertex_at(2)[2],
                                    after_io.at(i).vertex_at(2)[2], 1e-4);
    }
}
