#include <boost/mpl/list.hpp>

#include <boost/test/test_case_template.hpp>

#define BOOST_AUTO_TEST_MAIN
#include <boost/test/auto_unit_test.hpp>


#include "DynamicPriorityQueue.hpp"

typedef DynamicPriorityQueue<int>::ID ID;
typedef std::vector<ID> IDVector;


typedef DynamicPriorityQueue< int > IntegerDPQ;
typedef DynamicPriorityQueue< int, VolatileIDPolicy > VolatileIntegerDPQ;
typedef boost::mpl::list< IntegerDPQ, VolatileIntegerDPQ > both;
typedef boost::mpl::list< IntegerDPQ > novolatile;

BOOST_AUTO_TEST_CASE_TEMPLATE( testConstruction, DPQ, both )
{
    DPQ dpq;

    BOOST_CHECK( dpq.isEmpty() );
    BOOST_CHECK( dpq.checkConsistency() );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( testClear, DPQ, both )
{
    DPQ dpq;
    typedef typename DPQ::Index Index;

    dpq.push( 1 );
    dpq.push( 20 );
    dpq.push( 50 );

    BOOST_CHECK_EQUAL( Index( 3 ), dpq.getSize() );

    dpq.clear();

    BOOST_CHECK( dpq.isEmpty() );
    BOOST_CHECK( dpq.checkConsistency() );

    dpq.push( 2 );
    dpq.push( 20 );
    dpq.push( 30 );

    BOOST_CHECK_EQUAL( Index( 3 ), dpq.getSize() );

    dpq.clear();

    BOOST_CHECK( dpq.isEmpty() );
    BOOST_CHECK( dpq.checkConsistency() );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( testPush, DPQ, both )
{
    DPQ dpq;

    dpq.push( 1 );

    BOOST_CHECK( dpq.checkConsistency() );
    BOOST_CHECK( dpq.getTop() == 1.0 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testPushPop, DPQ, both )
{
    DynamicPriorityQueue<double> dpq;

    const ID id( dpq.push( 1 ) );

    BOOST_CHECK( dpq.checkConsistency() );
    BOOST_CHECK( dpq.getTop() == 1 );

    dpq.pop( id );

    BOOST_CHECK( dpq.checkConsistency() );
    BOOST_CHECK( dpq.isEmpty() );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testReplaceTop, DPQ, both )
{
    DPQ dpq;

    dpq.push( 4 );
    dpq.push( 2 );
    dpq.push( 1 );

    BOOST_CHECK_EQUAL( 1, dpq.getTop() );

    dpq.replaceTop( 3 );

    BOOST_CHECK( dpq.checkConsistency() );
    BOOST_CHECK_EQUAL( 2, dpq.getTop() );

    dpq.popTop();
    BOOST_CHECK_EQUAL( 3, dpq.getTop() );
    dpq.popTop();
    BOOST_CHECK_EQUAL( 4, dpq.getTop() );
    dpq.popTop();


    BOOST_CHECK( dpq.isEmpty() );
    BOOST_CHECK( dpq.checkConsistency() );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testReplace, DPQ, both )
{
    DPQ dpq;

    dpq.push( 5 );
    const ID id( dpq.push( 4 ) );
    dpq.push( 3 );
    dpq.push( 1 );

    BOOST_CHECK_EQUAL( 1, dpq.getTop() );

    dpq.replace( id, 2 );  // 4->2

    BOOST_CHECK( dpq.checkConsistency() );
    BOOST_CHECK_EQUAL( 1, dpq.getTop() );

    dpq.popTop();
    BOOST_CHECK_EQUAL( 2, dpq.getTop() );
    dpq.popTop();
    BOOST_CHECK_EQUAL( 3, dpq.getTop() );
    dpq.popTop();
    BOOST_CHECK_EQUAL( 5, dpq.getTop() );
    dpq.popTop();

    BOOST_CHECK( dpq.isEmpty() );
    BOOST_CHECK( dpq.checkConsistency() );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( testDuplicatedItems, DPQ, both )
{
    DPQ dpq;

    dpq.push( 1 );
    dpq.push( 2 );
    dpq.push( 1 );
    dpq.push( 2 );

    BOOST_CHECK( dpq.checkConsistency() );

    BOOST_CHECK( dpq.getTop() == 1 );
    dpq.popTop();
    BOOST_CHECK( dpq.getTop() == 1 );
    dpq.popTop();
    BOOST_CHECK( dpq.getTop() == 2 );
    dpq.popTop();
    BOOST_CHECK( dpq.getTop() == 2 );
    dpq.popTop();

    BOOST_CHECK( dpq.isEmpty() );

    BOOST_CHECK( dpq.checkConsistency() );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( testSimpleSorting, DPQ, both )
{
    DPQ dpq;

    const int MAXI( 100 );
    for( int i( MAXI ); i != 0  ; --i )
    {
        dpq.push( i );
    }

    BOOST_CHECK( dpq.checkConsistency() );

    int n( 0 );
    while( ! dpq.isEmpty() )
    {
        ++n;
        BOOST_CHECK_EQUAL( n, dpq.getTop() );
        dpq.popTop();
    }

    BOOST_CHECK_EQUAL( MAXI, n );

    BOOST_CHECK( dpq.isEmpty() );
    BOOST_CHECK( dpq.checkConsistency() );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( testSimpleSortingWithPops, DPQ, novolatile )
{
    DPQ dpq;
    typedef typename DPQ::Index Index;

    IDVector idVector;

    const Index MAXI( 100 );
    for( int n( MAXI ); n != 0  ; --n )
    {
        ID id( dpq.push( n ) );
        if( n == 11 || n == 45 )
        {
            idVector.push_back( id );
        }
    }

    BOOST_CHECK( dpq.checkConsistency() );

    BOOST_CHECK_EQUAL( MAXI, dpq.getSize() );

    for( IDVector::const_iterator i( idVector.begin() );
         i != idVector.end(); ++i )
    {
        dpq.pop( *i );
    }

    BOOST_CHECK_EQUAL( MAXI - 2, dpq.getSize() );

    int n( 0 );
    while( ! dpq.isEmpty() )
    {
        ++n;
        if( n == 11 || n == 45 )
        {
            continue; // skip
        }
        BOOST_CHECK_EQUAL( int( n ), dpq.getTop() );
        dpq.popTop();
    }

    BOOST_CHECK_EQUAL( MAXI, Index( n ) );

    BOOST_CHECK( dpq.isEmpty() );
    BOOST_CHECK( dpq.checkConsistency() );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( testInterleavedSorting, DPQ, both )
{
    DPQ dpq;
    typedef typename DPQ::Index Index;

    const Index MAXI( 101 );
    for( int i( MAXI-1 ); i != 0  ; i-=2 )
    {
        dpq.push( i );
    }

    for( int i( MAXI ); i != -1  ; i-=2 )
    {
        dpq.push( i );
    }

    BOOST_CHECK_EQUAL( MAXI, dpq.getSize() );

    BOOST_CHECK( dpq.checkConsistency() );

    int n( 0 );
    while( ! dpq.isEmpty() )
    {
        ++n;
        BOOST_CHECK_EQUAL( n, dpq.getTop() );
        dpq.popTop();
    }

    BOOST_CHECK_EQUAL( MAXI, Index( n ) );

    BOOST_CHECK( dpq.isEmpty() );
    BOOST_CHECK( dpq.checkConsistency() );
}


BOOST_AUTO_TEST_CASE_TEMPLATE( testInterleavedSortingWithPops, DPQ, 
                               novolatile )
{
    DPQ dpq;
    typedef typename DPQ::Index Index;

    IDVector idVector;

    const Index MAXI( 101 );
    for( int n( MAXI-1 ); n != 0  ; n-=2 )
    {
        const ID id( dpq.push( n ) );

        if( n == 12 || n == 46 )
        {
            idVector.push_back( id );
        }
    }

    dpq.pop( idVector.back() );
    idVector.pop_back();

    BOOST_CHECK_EQUAL( MAXI/2 -1, dpq.getSize() );

    BOOST_CHECK( dpq.checkConsistency() );

    for( int n( MAXI ); n != -1  ; n-=2 )
    {
        const ID id( dpq.push( n ) );

        if( n == 17 || n == 81 )
        {
            idVector.push_back( id );
        }
    }

    for( IDVector::const_iterator i( idVector.begin() );
         i != idVector.end(); ++i )
    {
        dpq.pop( *i );
    }

    BOOST_CHECK( dpq.checkConsistency() );
    BOOST_CHECK_EQUAL( MAXI-4, dpq.getSize() );

    int n( 0 );
    while( ! dpq.isEmpty() )
    {
        ++n;
        if( n == 12 || n == 46 || n == 17 || n == 81 )
        {
            continue;
        }
        BOOST_CHECK_EQUAL( n, dpq.getTop() );
        dpq.popTop();
    }

    BOOST_CHECK_EQUAL( MAXI, Index( n ) );

    BOOST_CHECK( dpq.isEmpty() );
    BOOST_CHECK( dpq.checkConsistency() );
}




