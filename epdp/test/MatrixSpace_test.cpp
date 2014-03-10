#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#define BOOST_TEST_MODULE "MatrixSpace_test"

#include <functional>
#include <iostream>
#include <cmath>
#include <set>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "Sphere.hpp"
#include "MatrixSpace.hpp"
#include "utils/random.hpp"
#include "GSLRandomNumberGenerator.hpp"

BOOST_AUTO_TEST_CASE(sized)
{
    typedef MatrixSpace<Sphere<double>, int> oc_type;
    typedef oc_type::position_type pos;
    BOOST_CHECK((is_sized<oc_type>::value));
    oc_type oc(1.0, 10);
    BOOST_CHECK_EQUAL(size(oc), 0);
    oc.update(std::make_pair(0, Sphere<double>(pos(0.2, 0.6, 0.4), 0.05)));
    BOOST_CHECK_EQUAL(size(oc), 1);
    oc.update(std::make_pair(1, Sphere<double>(pos(0.2, 0.6, 0.4), 0.05)));
    BOOST_CHECK_EQUAL(size(oc), 2);
    oc.update(std::make_pair(1, Sphere<double>(pos(0.1, 0.2, 0.3), 0.05)));
    BOOST_CHECK_EQUAL(size(oc), 2);
}

BOOST_AUTO_TEST_CASE(update)
{
    typedef MatrixSpace<Sphere<double>, int> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1.0, 10);
    BOOST_CHECK_CLOSE(0.1, oc.cell_size(), 0.001);

    {
        std::pair<oc_type::iterator, bool> ir(
                oc.update(std::make_pair(
                    0, oc_type::mapped_type(pos(0.2, 0.6, 0.4), 0.05))));
        BOOST_CHECK_EQUAL(true, ir.second);
        BOOST_CHECK(oc.end() != oc.find(0));
        BOOST_CHECK(oc.end() == oc.find(1));
    }
    {
        std::pair<oc_type::iterator, bool> ir(
                oc.update(std::make_pair(
                    0, oc_type::mapped_type(pos(0.2, 0.65, 0.4), 0.05))));
        BOOST_CHECK_EQUAL(false, ir.second);
        BOOST_CHECK_EQUAL(oc_type::mapped_type(pos(0.2, 0.65, 0.4), 0.05),
                (*ir.first).second);
        BOOST_CHECK(oc.end() != oc.find(0));
        BOOST_CHECK(oc.end() == oc.find(1));
    }
    {
        std::pair<oc_type::iterator, bool> ir(
                oc.update(std::make_pair(
                    0, oc_type::mapped_type(pos(0.2, 0.2, 0.4), 0.05))));
        BOOST_CHECK_EQUAL(false, ir.second);
        BOOST_CHECK_EQUAL(oc_type::mapped_type(pos(0.2, 0.2, 0.4), 0.05),
                (*ir.first).second);
        BOOST_CHECK(oc.end() != oc.find(0));
        BOOST_CHECK(oc.end() == oc.find(1));
    }
}

BOOST_AUTO_TEST_CASE(erase)
{
    typedef MatrixSpace<Sphere<double>, int> oc_type;
    typedef oc_type::position_type pos;

    for (int i = 0; i < 500; ++i)
    {
        oc_type oc(1.0, 10);

        if (i % 50 == 0)
        {
            std::cout << "*";
            std::cout.flush();
        }

        for (int j = 0; j < i; ++j)
        {
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(j));
            std::pair<oc_type::iterator, bool> ir(
                    oc.update(std::make_pair(
                        j, oc_type::mapped_type(pos(0.2 + 0.0001 * i, 0.6 + 0.0002 * i, 0.4), 0.001))));
            BOOST_CHECK_EQUAL(true, ir.second);
            for (int k = 0; k <= j; ++k)
                BOOST_CHECK(oc.end() != oc.find(k));
            BOOST_CHECK(oc.end() == oc.find(j + 1));
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(j + 1));
        }
        for (int j = i; --j >= 0;)
        {
            for (int k = 0; k <= j; ++k)
                BOOST_CHECK(oc.end() != oc.find(k));
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(j + 1));
            BOOST_CHECK(oc.erase(j));
            BOOST_CHECK(oc.end() == oc.find(j));
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(j));
        }
    }

    std::cout << std::endl;

    for (int i = 0; i < 500; ++i)
    {
        oc_type oc(1.0, 10);

        if (i % 50 == 0)
        {
            std::cout << "*";
            std::cout.flush();
        }

        for (int j = 0; j < i; ++j)
        {
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(j));
            std::pair<oc_type::iterator, bool> ir(
                    oc.update(std::make_pair(
                        j, oc_type::mapped_type(pos(0.2 + 0.0001 * i, 0.6 + 0.0002 * i, 0.4), 0.001))));
            BOOST_CHECK_EQUAL(true, ir.second);
            for (int k = 0; k <= j; ++k)
                BOOST_CHECK(oc.end() != oc.find(k));
            BOOST_CHECK(oc.end() == oc.find(j + 1));
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(j + 1));
        }
        for (int j = 0; j < i; ++j)
        {
            for (int k = j; k < i; ++k)
                BOOST_CHECK(oc.end() != oc.find(k));
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(i - j));
            BOOST_CHECK(oc.erase(j));
            BOOST_CHECK(oc.end() == oc.find(j));
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(i - j - 1));
        }
    }

    std::cout << std::endl;

    for (int i = 0; i < 500; ++i)
    {
        oc_type oc(1.0, 10);
        std::vector<int> id_list;

        if (i % 50 == 0)
        {
            std::cout << "*";
            std::cout.flush();
        }

        for (int j = 0; j < i; ++j)
        {
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(j));
            std::pair<oc_type::iterator, bool> ir(
                    oc.update(std::make_pair(
                        j, oc_type::mapped_type(pos(0.2 + 0.0001 * i, 0.6 + 0.0002 * i, 0.4), 0.001))));
            BOOST_CHECK_EQUAL(true, ir.second);
            for (int k = 0; k <= j; ++k)
                BOOST_CHECK(oc.end() != oc.find(k));
            BOOST_CHECK(oc.end() == oc.find(j + 1));
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(j + 1));
            id_list.push_back(j);
        }

        GSLRandomNumberGenerator rng;
        shuffle(rng, id_list);

        for (int j = 0; j < i; ++j)
        {
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(i - j));
            for (int k = 0; k < j; ++k)
                BOOST_CHECK(oc.end() == oc.find(id_list[k]));
            for (int k = j; k < i; ++k)
                BOOST_CHECK(oc.end() != oc.find(id_list[k]));
            BOOST_CHECK(oc.erase(id_list[j]));
            BOOST_CHECK(oc.size() == static_cast<oc_type::size_type>(i - j - 1));
        }
    }

    std::cout << std::endl;
}

template<typename Toc_>
struct collector
{
    collector(): result() {}

    void operator()(typename Toc_::iterator i,
            const typename Toc_::position_type&)
    {
        result.insert((*i).first);
    }

    std::set<typename Toc_::key_type> result;
};

template<typename Toc_>
struct collector2
{
    void operator()(typename Toc_::iterator i,
            const typename Toc_::position_type&)
    {
        result.insert((*i).first);
    }

    std::set<typename Toc_::key_type> result;
};

BOOST_AUTO_TEST_CASE(each_neighbor)
{
    typedef MatrixSpace<Sphere<double>, int> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1.0, 10);
    BOOST_CHECK_CLOSE(0.1, oc.cell_size(), 0.001);

    oc.update(std::make_pair(0, oc_type::mapped_type(pos(0.2, 0.6, 0.4), 0.05)));
    BOOST_CHECK(oc.end() != oc.find(0));
    BOOST_CHECK(oc.end() == oc.find(1));
    oc.update(std::make_pair(1, oc_type::mapped_type(pos(0.2, 0.7, 0.5), 0.05)));
    BOOST_CHECK(oc.end() != oc.find(0));
    BOOST_CHECK(oc.end() != oc.find(1));
    BOOST_CHECK(oc.end() == oc.find(2));
    oc.update(std::make_pair(2, oc_type::mapped_type(pos(0.9, 0.1, 0.4), 0.05)));
    BOOST_CHECK(oc.end() != oc.find(0));
    BOOST_CHECK(oc.end() != oc.find(1));
    BOOST_CHECK(oc.end() != oc.find(2));
    BOOST_CHECK(oc.end() == oc.find(3));
    oc.update(std::make_pair(3, oc_type::mapped_type(pos(0.9, 0.95, 0.4), 0.05)));
    BOOST_CHECK(oc.end() != oc.find(0));
    BOOST_CHECK(oc.end() != oc.find(1));
    BOOST_CHECK(oc.end() != oc.find(2));
    BOOST_CHECK(oc.end() != oc.find(3));
    BOOST_CHECK(oc.end() == oc.find(4));

    {
        collector<oc_type> col;
        oc.each_neighbor(oc.index(pos(0.2, 0.6, 0.4)), col);
        BOOST_CHECK_EQUAL(col.result.size(), 2);
        BOOST_CHECK(col.result.find(0) != col.result.end());
        BOOST_CHECK(col.result.find(1) != col.result.end());
    }

    {
        collector<oc_type> col;
        oc.each_neighbor(oc.index(pos(0.0, 0.1, 0.4)), col);
        BOOST_CHECK_EQUAL(col.result.size(), 0);
    }

    {
        collector2<oc_type> col2;
        oc.each_neighbor_cyclic(oc.index(pos(0.0, 0.1, 0.4)), col2);
        BOOST_CHECK_EQUAL(col2.result.size(), 1);
        BOOST_CHECK(col2.result.find(2) != col2.result.end());
    }
    {
        collector2<oc_type> col2;
        oc.each_neighbor_cyclic(oc.index(pos(0.9, 0.0, 0.4)), col2);
        BOOST_CHECK_EQUAL(col2.result.size(), 2);
        BOOST_CHECK(col2.result.find(2) != col2.result.end());
        BOOST_CHECK(col2.result.find(3) != col2.result.end());
    }
}

BOOST_AUTO_TEST_CASE(reupdate)
{
    typedef MatrixSpace<Sphere<double>, int> oc_type;
    typedef oc_type::position_type pos;
    oc_type oc(1.0, 10);
    BOOST_CHECK_CLOSE(0.1, oc.cell_size(), 0.001);

    BOOST_CHECK(oc.update(std::make_pair(0, oc_type::mapped_type(pos(0.2, 0.6, 0.4), 0.05))).second);
    BOOST_CHECK(oc.end() != oc.find(0));
    BOOST_CHECK(oc.end() == oc.find(1));

    {
        collector<oc_type> col;
        oc.each_neighbor(oc.index(pos(0.2, 0.6, 0.4)), col);
        BOOST_CHECK_EQUAL(col.result.size(), 1);
        BOOST_CHECK(col.result.find(0) != col.result.end());
        BOOST_CHECK(col.result.find(1) == col.result.end());
    }

    BOOST_CHECK(!oc.update(std::make_pair(0, oc_type::mapped_type(pos(0.13, 0.83, 0.43), 0.05))).second);
    BOOST_CHECK(oc.end() != oc.find(0));
    BOOST_CHECK(oc.end() == oc.find(1));

    {
        collector<oc_type> col;
        oc.each_neighbor(oc.index(pos(0.25, 0.62, 0.43)), col);
        BOOST_CHECK_EQUAL(col.result.size(), 0);
        BOOST_CHECK(col.result.find(0) == col.result.end());
        BOOST_CHECK(col.result.find(1) == col.result.end());
    }

    {
        collector<oc_type> col;
        oc.each_neighbor(oc.index(pos(0.23, 0.84, 0.45)), col);
        BOOST_CHECK_EQUAL(col.result.size(), 1);
        BOOST_CHECK(col.result.find(0) != col.result.end());
        BOOST_CHECK(col.result.find(1) == col.result.end());
    }
}

template<typename Toc_>
struct collector3
{
    collector3(): count(0) {}

    void operator()(typename Toc_::iterator i,
            typename Toc_::position_type const&)
    {
        ++count;
    }

    int count;
};

template<typename Toc_>
struct collector4
{
    collector4(): count(0) {}

    void operator()(typename Toc_::iterator i,
            typename Toc_::position_type const&)
    {
        ++count;
    }

    int count;
};

BOOST_AUTO_TEST_CASE(each_neighbor2)
{
    typedef MatrixSpace<Sphere<double>, int> oc_type;
    typedef oc_type::position_type pos;

    for (double r = 0.01; r < 0.1; r += 0.01)
    {
        std::cout << "*";
        std::cout.flush();
        for (double o = 0.0; o < 0.9; o += .001)
        {
            oc_type oc(1.0, 10);
            BOOST_CHECK_CLOSE(0.1, oc.cell_size(), 0.001);

            pos centre(o, o, o);

            for (int i = 0; i < 100; ++i)
            {
                double t1(M_PI * 2 * (i % 10) / 10.), t2(M_PI * 2 * (i / 10) / 10.);
                const double _x = cos(t1) * r;
                const pos p(centre[0] + _x * cos(t2),
                            centre[1] + sin(t1) * r,
                            centre[2] + _x * sin(t2));
                oc.update(std::make_pair(i, oc_type::mapped_type(p, r)));


                collector3<oc_type> col;
                oc.each_neighbor(oc.index(centre), col);
                BOOST_CHECK_EQUAL(col.count, i + 1);
            }
        }
    }
    std::cout << std::endl;

    for (double r = 0.01; r < 0.1; r += 0.01)
    {
        std::cout << "*";
        std::cout.flush();
        for (double o = 0.0; o < 0.9; o += .01)
        {
            oc_type oc(1.0, 10);
            BOOST_CHECK_CLOSE(0.1, oc.cell_size(), 0.001);

            pos centre(o, o, o);

            for (int i = 0; i < 100; ++i)
            {
                double t1(M_PI * 2 * (i % 10) / 10.), t2(M_PI * 2 * (i / 10) / 10.);
                const double _x = cos(t1) * r;
                const pos p(centre[0] + _x * cos(t2),
                            centre[1] + sin(t1) * r,
                            centre[2] + _x * sin(t2));
                oc.update(std::make_pair(i, oc_type::mapped_type(p, r)));


                collector4<oc_type> col;
                oc.each_neighbor_cyclic(oc.index(centre), col);
                BOOST_CHECK_EQUAL(col.count, i + 1);
            }
        }
    }
    std::cout << std::endl;
}

