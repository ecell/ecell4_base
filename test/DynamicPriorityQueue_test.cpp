
#define CPPUNIT_ENABLE_NAKED_ASSERT 1

#include <cppunit/TestFixture.h>
#include <cppunit/TestAssert.h>
#include <cppunit/extensions/HelperMacros.h>


#include "../DynamicPriorityQueue.hpp"

class DynamicPriorityQueueTest 
    : 
    public CppUnit::TestFixture 
{

public:

    typedef DynamicPriorityQueue< int > IntDPQ;
    typedef DynamicPriorityQueue< int, VolatileIDPolicy > VolatileIntDPQ;


    CPPUNIT_TEST_SUITE( DynamicPriorityQueueTest );

    CPPUNIT_TEST( testConstruction< IntDPQ > ); 
    CPPUNIT_TEST( testClear< IntDPQ > ); 
    CPPUNIT_TEST( testPush< IntDPQ > );
    CPPUNIT_TEST( testPushPop< IntDPQ > );
    CPPUNIT_TEST( testReplaceTop< IntDPQ > );
    CPPUNIT_TEST( testReplace< IntDPQ > );
    CPPUNIT_TEST( testDuplicatedItems< IntDPQ > );
    CPPUNIT_TEST( testSimpleSorting< IntDPQ > );
    CPPUNIT_TEST( testSimpleSortingWithPops< IntDPQ > );
    CPPUNIT_TEST( testInterleavedSorting< IntDPQ > );
    CPPUNIT_TEST( testInterleavedSortingWithPops< IntDPQ > );

    CPPUNIT_TEST( testConstruction< VolatileIntDPQ > ); 
    CPPUNIT_TEST( testClear< VolatileIntDPQ > ); 
    CPPUNIT_TEST( testPush< VolatileIntDPQ > );
    CPPUNIT_TEST( testPushPop< VolatileIntDPQ > );
    CPPUNIT_TEST( testReplaceTop< VolatileIntDPQ > );
    CPPUNIT_TEST( testReplace< VolatileIntDPQ > );
    CPPUNIT_TEST( testDuplicatedItems< VolatileIntDPQ > );
    CPPUNIT_TEST( testSimpleSorting< VolatileIntDPQ > );
    CPPUNIT_TEST( testInterleavedSorting< VolatileIntDPQ > );

    CPPUNIT_TEST_SUITE_END();

public:

    typedef DynamicPriorityQueue<int>::ID ID;
    typedef std::vector<ID> IDVector;

    void setUp()
    {
    }
    
    void tearDown() 
    {
    }


    template < class DPQ >
    void testConstruction()
    {
	DPQ dpq;

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

    template < class DPQ >
    void testClear()
    {
	DPQ dpq;

	dpq.push( 1 );
	dpq.push( 20 );
	dpq.push( 50 );

	//CPPUNIT_ASSERT_EQUAL( 3, dpq.getSize() );
	CPPUNIT_ASSERT( 3 == dpq.getSize() );

	dpq.clear();

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );

	dpq.push( 2 );
	dpq.push( 20 );
	dpq.push( 30 );

	CPPUNIT_ASSERT_EQUAL( 3, dpq.getSize() );

	dpq.clear();

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

    template < class DPQ >
    void testPush()
    {
	DPQ dpq;

	dpq.push( 1 );

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT( dpq.getTop() == 1.0 );
    }

    template < class DPQ >
    void testPushPop()
    {
	DynamicPriorityQueue<double> dpq;

	const ID id( dpq.push( 1 ) );

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT( dpq.getTop() == 1 );

	dpq.pop( id );

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT( dpq.isEmpty() );
    }

    template < class DPQ >
    void testReplaceTop()
    {
	DPQ dpq;

	dpq.push( 4 );
	dpq.push( 2 );
	dpq.push( 1 );

	CPPUNIT_ASSERT_EQUAL( 1, dpq.getTop() );

	dpq.replaceTop( 3 );

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT_EQUAL( 2, dpq.getTop() );

	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 3, dpq.getTop() );
	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 4, dpq.getTop() );
	dpq.popTop();


	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

    template < class DPQ >
    void testReplace()
    {
	DPQ dpq;

	dpq.push( 5 );
	const ID id( dpq.push( 4 ) );
	dpq.push( 3 );
	dpq.push( 1 );

	CPPUNIT_ASSERT_EQUAL( 1, dpq.getTop() );

	dpq.replace( id, 2 );  // 4->2

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT_EQUAL( 1, dpq.getTop() );

	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 2, dpq.getTop() );
	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 3, dpq.getTop() );
	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 5, dpq.getTop() );
	dpq.popTop();

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

    template < class DPQ >
    void testDuplicatedItems()
    {
	DPQ dpq;

	dpq.push( 1 );
	dpq.push( 2 );
	dpq.push( 1 );
	dpq.push( 2 );

	CPPUNIT_ASSERT( dpq.checkConsistency() );

	CPPUNIT_ASSERT( dpq.getTop() == 1 );
	dpq.popTop();
	CPPUNIT_ASSERT( dpq.getTop() == 1 );
	dpq.popTop();
	CPPUNIT_ASSERT( dpq.getTop() == 2 );
	dpq.popTop();
	CPPUNIT_ASSERT( dpq.getTop() == 2 );
	dpq.popTop();

	CPPUNIT_ASSERT( dpq.isEmpty() );

	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }


    template < class DPQ >
    void testSimpleSorting()
    {
	DPQ dpq;

	const int MAXI( 100 );
	for( int i( MAXI ); i != 0  ; --i )
	{
	    dpq.push( i );
	}

	CPPUNIT_ASSERT( dpq.checkConsistency() );

	int n( 0 );
	while( ! dpq.isEmpty() )
	{
	    ++n;
	    CPPUNIT_ASSERT_EQUAL( n, dpq.getTop() );
	    dpq.popTop();
	}

	CPPUNIT_ASSERT_EQUAL( MAXI, n );

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }


    template < class DPQ >
    void testSimpleSortingWithPops()
    {
	DPQ dpq;

	IDVector idVector;

	const int MAXI( 100 );
	for( int n( MAXI ); n != 0  ; --n )
	{
	    ID id( dpq.push( n ) );
	    if( n == 11 || n == 45 )
	    {
		idVector.push_back( id );
	    }
	}

	CPPUNIT_ASSERT( dpq.checkConsistency() );

	CPPUNIT_ASSERT_EQUAL( MAXI, dpq.getSize() );

	for( IDVector::const_iterator i( idVector.begin() );
	     i != idVector.end(); ++i )
	{
	    dpq.pop( *i );
	}

	CPPUNIT_ASSERT_EQUAL( MAXI - 2, dpq.getSize() );

	int n( 0 );
	while( ! dpq.isEmpty() )
	{
	    ++n;
	    if( n == 11 || n == 45 )
	    {
		continue; // skip
	    }
	    CPPUNIT_ASSERT_EQUAL( n, dpq.getTop() );
	    dpq.popTop();
	}

	CPPUNIT_ASSERT_EQUAL( MAXI, n );

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }


    template < class DPQ >
    void testInterleavedSorting()
    {
	DPQ dpq;

	const int MAXI( 101 );
	for( int i( MAXI-1 ); i != 0  ; i-=2 )
	{
	    dpq.push( i );
	}

	for( int i( MAXI ); i != -1  ; i-=2 )
	{
	    dpq.push( i );
	}

	CPPUNIT_ASSERT_EQUAL( MAXI, dpq.getSize() );

	CPPUNIT_ASSERT( dpq.checkConsistency() );

	int n( 0 );
	while( ! dpq.isEmpty() )
	{
	    ++n;
	    CPPUNIT_ASSERT_EQUAL( n, dpq.getTop() );
	    dpq.popTop();
	}

	CPPUNIT_ASSERT_EQUAL( MAXI, n );

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }


    template < class DPQ >
    void testInterleavedSortingWithPops()
    {
	DPQ dpq;

	IDVector idVector;

	const int MAXI( 101 );
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

	CPPUNIT_ASSERT_EQUAL( MAXI/2 -1, dpq.getSize() );

	CPPUNIT_ASSERT( dpq.checkConsistency() );

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

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT_EQUAL( MAXI-4, dpq.getSize() );

	int n( 0 );
	while( ! dpq.isEmpty() )
	{
	    ++n;
	    if( n == 12 || n == 46 || n == 17 || n == 81 )
	    {
		continue;
	    }
	    CPPUNIT_ASSERT_EQUAL( n, dpq.getTop() );
	    dpq.popTop();
	}

	CPPUNIT_ASSERT_EQUAL( MAXI, n );

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( DynamicPriorityQueueTest );

