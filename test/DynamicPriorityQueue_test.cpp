#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>


#include "../DynamicPriorityQueue.hpp"

class DynamicPriorityQueueTest 
    : 
    public CppUnit::TestFixture 
{

public:

    CPPUNIT_TEST_SUITE( DynamicPriorityQueueTest );
    CPPUNIT_TEST( testConstruction ); 
    CPPUNIT_TEST( testPush );
    CPPUNIT_TEST( testPushPop );
    CPPUNIT_TEST( testReplaceTop );
    CPPUNIT_TEST( testReplace );
    CPPUNIT_TEST( testDuplicatedItems );
    CPPUNIT_TEST( testSimpleSorting );
    CPPUNIT_TEST( testSimpleSortingWithPops );
    CPPUNIT_TEST( testInterleavedSorting );
    CPPUNIT_TEST( testInterleavedSortingWithPops );
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
    
    void testConstruction()
    {
	DynamicPriorityQueue<double> dpq;

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

    void testPush()
    {
	DynamicPriorityQueue<double> dpq;

	dpq.pushItem( 1.0 );

	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT( dpq.getTopItem() == 1.0 );
    }

    void testPushPop()
    {
	DynamicPriorityQueue<double> dpq;

	const ID id( dpq.pushItem( 1.0 ) );

	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT( dpq.getTopItem() == 1.0 );

	dpq.popItem( id );

	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT( dpq.isEmpty() );
    }

    void testReplaceTop()
    {
	DynamicPriorityQueue<int> dpq;

	dpq.pushItem( 4 );
	dpq.pushItem( 2 );
	dpq.pushItem( 1 );

	CPPUNIT_ASSERT_EQUAL( 1, dpq.getTopItem() );

	dpq.replaceTop( 3 );

	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
	CPPUNIT_ASSERT_EQUAL( 2, dpq.getTopItem() );

	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 3, dpq.getTopItem() );
	dpq.popTop();
	CPPUNIT_ASSERT_EQUAL( 4, dpq.getTopItem() );
	dpq.popTop();


	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

    void testReplace()
    {
	DynamicPriorityQueue<int> dpq;

	dpq.pushItem( 5 );
	const ID id( dpq.pushItem( 4 ) );
	dpq.pushItem( 3 );
	dpq.pushItem( 1 );

	CPPUNIT_ASSERT_EQUAL( 1, dpq.getTopItem() );

	dpq.replaceItem( id, 2 );  // 4->2

	CPPUNIT_ASSERT( dpq.checkSize() );
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
	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

    void testDuplicatedItems()
    {
	DynamicPriorityQueue<int> dpq;

	dpq.pushItem( 1 );
	dpq.pushItem( 2 );
	dpq.pushItem( 1 );
	dpq.pushItem( 2 );

	CPPUNIT_ASSERT( dpq.checkSize() );
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

	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }


    void testSimpleSorting()
    {
	DynamicPriorityQueue<int> dpq;

	const int MAXI( 100 );
	for( int i( MAXI ); i != 0  ; --i )
	{
	    dpq.pushItem( i );
	}

	CPPUNIT_ASSERT( dpq.checkSize() );
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
	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }


    void testSimpleSortingWithPops()
    {
	DynamicPriorityQueue<int> dpq;

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

	CPPUNIT_ASSERT( dpq.checkSize() );
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
	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }


    void testInterleavedSorting()
    {
	DynamicPriorityQueue<int> dpq;

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

	CPPUNIT_ASSERT( dpq.checkSize() );
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
	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }


    void testInterleavedSortingWithPops()
    {
	DynamicPriorityQueue<int> dpq;

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

	CPPUNIT_ASSERT( dpq.checkSize() );
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

	CPPUNIT_ASSERT( dpq.checkSize() );
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
	CPPUNIT_ASSERT( dpq.checkSize() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( DynamicPriorityQueueTest );

