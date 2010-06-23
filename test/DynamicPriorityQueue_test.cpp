#define BOOST_TEST_MODULE "DynamicPriorityQueue"

#include <boost/mpl/list.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/test/test_case_template.hpp>

#include "DynamicPriorityQueue.hpp"

typedef DynamicPriorityQueue<int>::identifier_type identifier_type;
typedef std::vector<identifier_type> identifier_vector;


typedef DynamicPriorityQueue<int > IntegerDPQ;
typedef DynamicPriorityQueue<int, std::less_equal<int>, volatile_id_policy<> > VolatileIntegerDPQ;
typedef boost::mpl::list<IntegerDPQ, VolatileIntegerDPQ> both;
typedef boost::mpl::list<IntegerDPQ> novolatile;

BOOST_AUTO_TEST_CASE_TEMPLATE(testConstruction, DPQ, both)
{
    DPQ dpq;

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testClear, DPQ, both)
{
    DPQ dpq;
    typedef typename DPQ::index_type Index;

    dpq.push(1);
    dpq.push(20);
    dpq.push(50);

    BOOST_CHECK_EQUAL(Index(3), dpq.size());

    dpq.clear();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());

    dpq.push(2);
    dpq.push(20);
    dpq.push(30);

    BOOST_CHECK_EQUAL(Index(3), dpq.size());

    dpq.clear();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testClear_manyItems, DPQ, both)
{
    DPQ dpq;
    typedef typename DPQ::index_type Index;

    for (int i(0); i < 70000; ++i)
    {
        dpq.push(i);
    }

    BOOST_CHECK(dpq.check());

    dpq.clear();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());

    dpq.push(2);
    dpq.push(20);
    dpq.push(30);

    BOOST_CHECK_EQUAL(Index(3), dpq.size());

    dpq.clear();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}



BOOST_AUTO_TEST_CASE_TEMPLATE(testPush, DPQ, both)
{
    DPQ dpq;

    dpq.push(1);

    BOOST_CHECK(dpq.check());
    BOOST_CHECK(dpq.top().second == 1.0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(testPushPop, DPQ, both)
{
    DynamicPriorityQueue<double> dpq;

    const identifier_type id(dpq.push(1));

    BOOST_CHECK(dpq.check());
    BOOST_CHECK(dpq.top().second == 1);

    dpq.pop(id);

    BOOST_CHECK(dpq.check());
    BOOST_CHECK(dpq.empty());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testSecond01, DPQ, both)
{
    DPQ dpq;

    BOOST_CHECK_THROW(dpq.second(), std::out_of_range);

    dpq.push(4);

    BOOST_CHECK_EQUAL(4, dpq.top().second);
    BOOST_CHECK_THROW(dpq.second(), std::out_of_range);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(testSecond2, DPQ, both)
{
    {
        DPQ dpq;

        dpq.push(1);
        dpq.push(4);

        BOOST_CHECK_EQUAL(1, dpq.top().second);
        BOOST_CHECK_EQUAL(4, dpq.second().second);
    }


    {
        DPQ dpq;

        dpq.push(4);
        dpq.push(1);

        BOOST_CHECK_EQUAL(1, dpq.top().second);
        BOOST_CHECK_EQUAL(4, dpq.second().second);
    }
}



BOOST_AUTO_TEST_CASE_TEMPLATE(testSecondN, DPQ, both)
{
    DPQ dpq;

    dpq.push(2);
    dpq.push(4);
    dpq.push(1);
    dpq.push(5);


    BOOST_CHECK_EQUAL(1, dpq.top().second);
    BOOST_CHECK_EQUAL(2, dpq.second().second);

    dpq.replace(std::make_pair(dpq.top().first, 3));

    BOOST_CHECK(dpq.check());
    BOOST_CHECK_EQUAL(2, dpq.top().second);
    BOOST_CHECK_EQUAL(3, dpq.second().second);

    dpq.pop();
    BOOST_CHECK_EQUAL(3, dpq.top().second);
    BOOST_CHECK_EQUAL(4, dpq.second().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(4, dpq.top().second);
    BOOST_CHECK_EQUAL(5, dpq.second().second);
    dpq.pop();

    BOOST_CHECK_EQUAL(5, dpq.top().second);
    dpq.pop();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(testReplace1, DPQ, both)
{
    DPQ dpq;

    dpq.push(5);      
    const identifier_type id(dpq.push(4));     
    dpq.push(3);      
    dpq.push(1);      
    BOOST_CHECK(dpq.check());

    BOOST_CHECK_EQUAL(4, dpq.size());
    BOOST_CHECK_EQUAL(1, dpq.top().second);

    dpq.replace(std::make_pair(id, 2));  // 4->2 up

    BOOST_CHECK(dpq.check());
    BOOST_CHECK_EQUAL(4, dpq.size());

    BOOST_CHECK_EQUAL(1, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(2, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(3, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(5, dpq.top().second);
    dpq.pop();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(testReplace2, DPQ, both)
{
    DPQ dpq;

    dpq.push(5);
    dpq.push(4);
    const identifier_type id(dpq.push(3));
    dpq.push(1);
    BOOST_CHECK(dpq.check());

    BOOST_CHECK_EQUAL(4, dpq.size());
    BOOST_CHECK_EQUAL(1, dpq.top().second);

    dpq.replace(std::make_pair(id, 6));  // 3->6 down

    BOOST_CHECK(dpq.check());
    BOOST_CHECK_EQUAL(4, dpq.size());

    BOOST_CHECK_EQUAL(1, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(4, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(5, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(6, dpq.top().second);
    dpq.pop();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(testReplace3, DPQ, both)
{
    DPQ dpq;

    dpq.push(5);
    dpq.push(4);
    const identifier_type id(dpq.push(3));
    dpq.push(1);
    BOOST_CHECK(dpq.check());

    BOOST_CHECK_EQUAL(4, dpq.size());
    BOOST_CHECK_EQUAL(1, dpq.top().second);

    dpq.replace(std::make_pair(id, 3));  // 3->3 no change

    BOOST_CHECK(dpq.check());
    BOOST_CHECK_EQUAL(4, dpq.size());

    BOOST_CHECK_EQUAL(1, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(3, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(4, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(5, dpq.top().second);
    dpq.pop();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testReplace4, DPQ, both)
{
    DPQ dpq;

    dpq.push(5);
    dpq.push(4);
    dpq.push(3);
    const identifier_type id(dpq.push(1));
    BOOST_CHECK(dpq.check());

    BOOST_CHECK_EQUAL(4, dpq.size());
    BOOST_CHECK_EQUAL(1, dpq.top().second);

    dpq.replace(std::make_pair(id, 1));  // 1->1 top no change

    BOOST_CHECK(dpq.check());
    BOOST_CHECK_EQUAL(4, dpq.size());

    BOOST_CHECK_EQUAL(1, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(3, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(4, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(5, dpq.top().second);
    dpq.pop();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testReplace5, DPQ, both)
{
    DPQ dpq;

    dpq.push(5);
    dpq.push(4);
    dpq.push(3);
    const identifier_type id(dpq.push(1));
    BOOST_CHECK(dpq.check());

    BOOST_CHECK_EQUAL(4, dpq.size());
    BOOST_CHECK_EQUAL(1, dpq.top().second);

    dpq.replace(std::make_pair(id, 0));  // 1->0 top up, no change

    BOOST_CHECK(dpq.check());
    BOOST_CHECK_EQUAL(4, dpq.size());

    BOOST_CHECK_EQUAL(0, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(3, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(4, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(5, dpq.top().second);
    dpq.pop();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(testReplace6, DPQ, both)
{
    DPQ dpq;

    dpq.push(5);
    dpq.push(4);
    dpq.push(3);
    const identifier_type id(dpq.push(1));
    BOOST_CHECK(dpq.check());

    BOOST_CHECK_EQUAL(4, dpq.size());
    BOOST_CHECK_EQUAL(1, dpq.top().second);

    dpq.replace(std::make_pair(id, 3));  // 1->3 top down

    BOOST_CHECK(dpq.check());
    BOOST_CHECK_EQUAL(4, dpq.size());

    BOOST_CHECK_EQUAL(3, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(3, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(4, dpq.top().second);
    dpq.pop();
    BOOST_CHECK_EQUAL(5, dpq.top().second);
    dpq.pop();

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testDuplicatedItems, DPQ, both)
{
    DPQ dpq;

    dpq.push(1);
    dpq.push(2);
    dpq.push(1);
    dpq.push(2);

    BOOST_CHECK(dpq.check());

    BOOST_CHECK(dpq.top().second == 1);
    dpq.pop();
    BOOST_CHECK(dpq.top().second == 1);
    dpq.pop();
    BOOST_CHECK(dpq.top().second == 2);
    dpq.pop();
    BOOST_CHECK(dpq.top().second == 2);
    dpq.pop();

    BOOST_CHECK(dpq.empty());

    BOOST_CHECK(dpq.check());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testPushAfterPop, DPQ, both)
{
    DPQ dpq;

    dpq.push(1);
    const identifier_type id1(dpq.push(2));
    dpq.push(3);
    dpq.push(4);
    BOOST_CHECK(dpq.size() == 4);

    BOOST_CHECK(dpq.check());

    BOOST_CHECK(dpq.top().second == 1);
    dpq.pop(id1);
    BOOST_CHECK(dpq.size() == 3);
    BOOST_CHECK(dpq.check());

    dpq.push(1);
    BOOST_CHECK(dpq.size() == 4);
    BOOST_CHECK(dpq.check());

    BOOST_CHECK(dpq.top().second == 1);
    dpq.pop();
    BOOST_CHECK(dpq.top().second == 1);
    dpq.pop();
    BOOST_CHECK(dpq.top().second == 3);
    dpq.pop();
    BOOST_CHECK(dpq.top().second == 4);
    dpq.pop();
    BOOST_CHECK(dpq.empty());

    BOOST_CHECK(dpq.check());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testPushAfterPopTwoItems, DPQ, both)
{
    DPQ dpq;

    const identifier_type id1(dpq.push(2));
    dpq.push(1);
    BOOST_CHECK(dpq.size() == 2);

    BOOST_CHECK(dpq.check());

    BOOST_CHECK(dpq.top().second == 1);
    dpq.pop(id1);
    BOOST_CHECK(dpq.size() == 1);
    BOOST_CHECK(dpq.check());
}



BOOST_AUTO_TEST_CASE_TEMPLATE(testSimpleSorting, DPQ, both)
{
    DPQ dpq;

    const int MAXI(100);
    for (int i(MAXI); i != 0  ; --i)
    {
        dpq.push(i);
    }

    BOOST_CHECK(dpq.check());

    int n(0);
    while (! dpq.empty())
    {
        ++n;
        BOOST_CHECK_EQUAL(n, dpq.top().second);
        dpq.pop();
    }

    BOOST_CHECK_EQUAL(MAXI, n);

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testSimpleSortingWithPops, DPQ, novolatile)
{
    DPQ dpq;
    typedef typename DPQ::index_type Index;

    identifier_vector idVector;

    const Index MAXI(100);
    for (int n(MAXI); n != 0  ; --n)
    {
        identifier_type id(dpq.push(n));
        if (n == 11 || n == 45)
        {
            idVector.push_back(id);
        }
    }

    BOOST_CHECK(dpq.check());

    BOOST_CHECK_EQUAL(MAXI, dpq.size());

    for (identifier_vector::const_iterator i(idVector.begin());
         i != idVector.end(); ++i)
    {
        dpq.pop(*i);
    }

    BOOST_CHECK_EQUAL(MAXI - 2, dpq.size());

    int n(0);
    while (! dpq.empty())
    {
        ++n;
        if (n == 11 || n == 45)
        {
            continue; // skip
        }
        BOOST_CHECK_EQUAL(int(n), dpq.top().second);
        dpq.pop();
    }

    BOOST_CHECK_EQUAL(MAXI, Index(n));

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testInterleavedSorting, DPQ, both)
{
    DPQ dpq;
    typedef typename DPQ::index_type Index;

    const Index MAXI(101);
    for (int i(MAXI-1); i != 0  ; i-=2)
    {
        dpq.push(i);
    }

    for (int i(MAXI); i != -1  ; i-=2)
    {
        dpq.push(i);
    }

    BOOST_CHECK_EQUAL(MAXI, dpq.size());

    BOOST_CHECK(dpq.check());

    int n(0);
    while (! dpq.empty())
    {
        ++n;
        BOOST_CHECK_EQUAL(n, dpq.top().second);
        dpq.pop();
    }

    BOOST_CHECK_EQUAL(MAXI, Index(n));

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}


BOOST_AUTO_TEST_CASE_TEMPLATE(testInterleavedSortingWithPops, DPQ, 
                               novolatile)
{
    DPQ dpq;
    typedef typename DPQ::index_type Index;

    identifier_vector idVector;

    const Index MAXI(101);
    for (int n(MAXI-1); n != 0  ; n-=2)
    {
        const identifier_type id(dpq.push(n));

        if (n == 12 || n == 46)
        {
            idVector.push_back(id);
        }
    }

    dpq.pop(idVector.back());
    idVector.pop_back();

    BOOST_CHECK_EQUAL(MAXI/2 -1, dpq.size());

    BOOST_CHECK(dpq.check());

    for (int n(MAXI); n != -1  ; n-=2)
    {
        const identifier_type id(dpq.push(n));

        if (n == 17 || n == 81)
        {
            idVector.push_back(id);
        }
    }

    for (identifier_vector::const_iterator i(idVector.begin());
         i != idVector.end(); ++i)
    {
        dpq.pop(*i);
    }

    BOOST_CHECK(dpq.check());
    BOOST_CHECK_EQUAL(MAXI-4, dpq.size());

    int n(0);
    while (! dpq.empty())
    {
        ++n;
        if (n == 12 || n == 46 || n == 17 || n == 81)
        {
            continue;
        }
        BOOST_CHECK_EQUAL(n, dpq.top().second);
        dpq.pop();
    }

    BOOST_CHECK_EQUAL(MAXI, Index(n));

    BOOST_CHECK(dpq.empty());
    BOOST_CHECK(dpq.check());
}
