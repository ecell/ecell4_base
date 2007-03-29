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

    typedef DynamicPriorityQueue< int > IntegerDPQ;
    typedef DynamicPriorityQueue< int, VolatileIDPolicy > VolatileIntegerDPQ;


    CPPUNIT_TEST_SUITE( DynamicPriorityQueueTest );

    CPPUNIT_TEST( testConstruction< IntegerDPQ > ); 
    CPPUNIT_TEST( testClear< IntegerDPQ > ); 
    CPPUNIT_TEST( testPush< IntegerDPQ > );
    CPPUNIT_TEST( testPushPop< IntegerDPQ > );
    CPPUNIT_TEST( testReplaceTop< IntegerDPQ > );
    CPPUNIT_TEST( testReplace< IntegerDPQ > );
    CPPUNIT_TEST( testDuplicatedItems< IntegerDPQ > );
    CPPUNIT_TEST( testSimpleSorting< IntegerDPQ > );
    CPPUNIT_TEST( testSimpleSortingWithPops< IntegerDPQ > );
    CPPUNIT_TEST( testIntegererleavedSorting< IntegerDPQ > );
    CPPUNIT_TEST( testIntegererleavedSortingWithPops< IntegerDPQ > );

    CPPUNIT_TEST( testConstruction< VolatileIntegerDPQ > ); 
    CPPUNIT_TEST( testClear< VolatileIntegerDPQ > ); 
    CPPUNIT_TEST( testPush< VolatileIntegerDPQ > );
    CPPUNIT_TEST( testPushPop< VolatileIntegerDPQ > );
    CPPUNIT_TEST( testReplaceTop< VolatileIntegerDPQ > );
    CPPUNIT_TEST( testReplace< VolatileIntegerDPQ > );
    CPPUNIT_TEST( testDuplicatedItems< VolatileIntegerDPQ > );
    CPPUNIT_TEST( testSimpleSorting< VolatileIntegerDPQ > );
    CPPUNIT_TEST( testIntegererleavedSorting< VolatileIntegerDPQ > );

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
	typedef typename DPQ::Index Index;

	dpq.push( 1 );
	dpq.push( 20 );
	dpq.push( 50 );

	CPPUNIT_ASSERT_EQUAL( Index( 3 ), dpq.getSize() );

	dpq.clear();

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );

	dpq.push( 2 );
	dpq.push( 20 );
	dpq.push( 30 );

	CPPUNIT_ASSERT_EQUAL( Index( 3 ), dpq.getSize() );

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
	    CPPUNIT_ASSERT_EQUAL( int( n ), dpq.getTop() );
	    dpq.popTop();
	}

	CPPUNIT_ASSERT_EQUAL( MAXI, Index( n ) );

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }


    template < class DPQ >
    void testIntegererleavedSorting()
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

	CPPUNIT_ASSERT_EQUAL( MAXI, dpq.getSize() );

	CPPUNIT_ASSERT( dpq.checkConsistency() );

	int n( 0 );
	while( ! dpq.isEmpty() )
	{
	    ++n;
	    CPPUNIT_ASSERT_EQUAL( n, dpq.getTop() );
	    dpq.popTop();
	}

	CPPUNIT_ASSERT_EQUAL( MAXI, Index( n ) );

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }


    template < class DPQ >
    void testIntegererleavedSortingWithPops()
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

	CPPUNIT_ASSERT_EQUAL( MAXI, Index( n ) );

	CPPUNIT_ASSERT( dpq.isEmpty() );
	CPPUNIT_ASSERT( dpq.checkConsistency() );
    }

};

CPPUNIT_TEST_SUITE_REGISTRATION( DynamicPriorityQueueTest );

