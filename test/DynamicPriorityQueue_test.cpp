#include <cppunit/TestFixture.h>
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

	dpq.pushItem( 1 );
	dpq.pushItem( 20 );
	dpq.pushItem( 50 );

	CPPUNIT_ASSERT_EQUAL( 3, dpq.getSize() );

	dpq.clear();

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );

	dpq.pushItem( 2 );
	dpq.pushItem( 20 );
	dpq.pushItem( 30 );

	CPPUNIT_ASSERT_EQUAL( 3, dpq.getSize() );

	dpq.clear();

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

    template < class DPQ >
    void testPush()
    {
	DPQ dpq;

	dpq.pushItem( 1 );

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT( dpq.getTopItem() == 1.0 );
    }

    template < class DPQ >
    void testPushPop()
    {
	DynamicPriorityQueue<double> dpq;

	const ID id( dpq.pushItem( 1 ) );

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT( dpq.getTopItem() == 1 );

	dpq.popItem( id );

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT( dpq.isEmpty() );
    }

    template < class DPQ >
    void testReplaceTop()
    {
	DPQ dpq;

	dpq.pushItem( 4 );
	dpq.pushItem( 2 );
	dpq.pushItem( 1 );

	CPPUNIT_ASSERT_EQUAL( 1, dpq.getTopItem() );

	dpq.replaceTop( 3 );

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT_EQUAL( 2, dpq.getTopItem() );

	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 3, dpq.getTopItem() );
	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 4, dpq.getTopItem() );
	dpq.popTop();


	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

    template < class DPQ >
    void testReplace()
    {
	DPQ dpq;

	dpq.pushItem( 5 );
	const ID id( dpq.pushItem( 4 ) );
	dpq.pushItem( 3 );
	dpq.pushItem( 1 );

	CPPUNIT_ASSERT_EQUAL( 1, dpq.getTopItem() );

	dpq.replaceItem( id, 2 );  // 4->2

	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT_EQUAL( 1, dpq.getTopItem() );

	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 2, dpq.getTopItem() );
	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 3, dpq.getTopItem() );
	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 5, dpq.getTopItem() );
	dpq.popTop();

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

    template < class DPQ >
    void testDuplicatedItems()
    {
	DPQ dpq;

	dpq.pushItem( 1 );
	dpq.pushItem( 2 );
	dpq.pushItem( 1 );
	dpq.pushItem( 2 );

	CPPUNIT_ASSERT( dpq.checkConsistency() );

	CPPUNIT_ASSERT( dpq.getTopItem() == 1 );
	dpq.popTop();
	CPPUNIT_ASSERT( dpq.getTopItem() == 1 );
	dpq.popTop();
	CPPUNIT_ASSERT( dpq.getTopItem() == 2 );
	dpq.popTop();
	CPPUNIT_ASSERT( dpq.getTopItem() == 2 );
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
	    dpq.pushItem( i );
	}

	CPPUNIT_ASSERT( dpq.checkConsistency() );

	int n( 0 );
	while( ! dpq.isEmpty() )
	{
	    ++n;
	    CPPUNIT_ASSERT_EQUAL( n, dpq.getTopItem() );
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
	    ID id( dpq.pushItem( n ) );
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
	    dpq.popItem( *i );
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
	    CPPUNIT_ASSERT_EQUAL( n, dpq.getTopItem() );
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
	    dpq.pushItem( i );
	}

	for( int i( MAXI ); i != -1  ; i-=2 )
	{
	    dpq.pushItem( i );
	}

	CPPUNIT_ASSERT_EQUAL( MAXI, dpq.getSize() );

	CPPUNIT_ASSERT( dpq.checkConsistency() );

	int n( 0 );
	while( ! dpq.isEmpty() )
	{
	    ++n;
	    CPPUNIT_ASSERT_EQUAL( n, dpq.getTopItem() );
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
	    const ID id( dpq.pushItem( n ) );

	    if( n == 12 || n == 46 )
	    {
		idVector.push_back( id );
	    }
	}

	dpq.popItem( idVector.back() );
	idVector.pop_back();

	CPPUNIT_ASSERT_EQUAL( MAXI/2 -1, dpq.getSize() );

	CPPUNIT_ASSERT( dpq.checkConsistency() );

	for( int n( MAXI ); n != -1  ; n-=2 )
	{
	    const ID id( dpq.pushItem( n ) );

	    if( n == 17 || n == 81 )
	    {
		idVector.push_back( id );
	    }
	}

	for( IDVector::const_iterator i( idVector.begin() );
	     i != idVector.end(); ++i )
	{
	    dpq.popItem( *i );
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
	    CPPUNIT_ASSERT_EQUAL( n, dpq.getTopItem() );
	    dpq.popTop();
	}

	CPPUNIT_ASSERT_EQUAL( MAXI, n );

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( DynamicPriorityQueueTest );

