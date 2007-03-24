//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 1996-2005 Keio University
//                Copyright (C) 2005-2007 The Molecular Sciences Institute
//
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//
// E-Cell is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
// 
// E-Cell is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public
// License along with E-Cell -- see the file COPYING.
// If not, write to the Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
//END_HEADER
//
// written by Eiichiro Adachi and Koichi Takahashi
//

//#define NDEBUG
#undef NDEBUG

#ifndef __DYNAMICPRIORITYQUEUE_HPP
#define __DYNAMICPRIORITYQUEUE_HPP
#include <assert.h>
#include <vector>
#include <algorithm>
#include <map>
#include <tr1/unordered_map>


template < class T >
struct PtrGreater
{
    bool operator()( T x, T y ) const { return *y < *x; }
};


template < typename Item >
class DynamicPriorityQueue
{
  

public:

    typedef long long unsigned int ID;

    typedef std::vector< Item >    ItemVector;
    //typedef std::vector< Item* >   ItemPtrVector;

    typedef typename ItemVector::difference_type       Index;

    typedef std::vector< Index >  IndexVector;
    typedef std::vector< ID >     IDVector;

//    typedef std::tr1::unordered_map<const ID, Index> IndexMap;
    typedef std::map<const ID, Index> IndexMap;


    DynamicPriorityQueue();
  
    void move( const Index index )
    {
	const Index pos( this->positionVector[index] );
	movePos( pos );
    }

    inline void movePos( const Index pos );

    void moveTop()
    {
	moveDownPos( 0 );
    }

    const Index getIndex( const ID id ) const
    {
	return this->indexMap.at( id );
    }

    const Index getIDByIndex( const Index index ) const
    {
	return this->idVector[ index ];
    }

    const Index getTopIndex() const 
    {
	assert( positionVector[this->heap[0]] == 0 );
	assert( this->heap[0] <= getSize() );
	return this->heap[0];
    }

    const Item& getTopItem() const
    {
	return this->itemVector[ this->heap[0] ];
    }

    Item& getTopItem()
    {
	return this->itemVector[ this->heap[0] ];
    }

    const Item& getItemByIndex( const Index index ) const
    {
	return this->itemVector[ index ];
    }

    Item& getItemByIndex( const Index index )
    {
	return this->itemVector[ index ];
    }

    const Item& getItem( const ID id ) const
    {
	return this->getItemByIndex( this->indexMap[ id ] );
    }

    Item& getItem( const ID id )
    {
	return this->getItemByIndex( this->indexMap[ id ] );
    }


    void popTop()
    {
	return popItemByIndex( getTopIndex() );
    }

    void popItem( const ID id )
    {
	return popItemByIndex( getIndex( id ) );
    }

    void popItemByIndex( const Index index );

    const ID pushItem( const Item& item )
    {
	consistencyCheck();

	const Index index( this->size );
    
	++this->size;
    
//	this->itemVector.push_back( item );
	// index == pos at this time.
//	this->heap.push_back( index );
//	this->positionVector.push_back( index );

	this->itemVector.resize( getSize() );
	this->heap.resize( getSize() );
	this->positionVector.resize( getSize() );
	this->idVector.resize( getSize() );

	this->itemVector[index] = item;
	// index == pos at this time.
	this->heap[index] = index;
	this->positionVector[index] = index;

	const ID id( this->idCounter );
	++this->idCounter;

	this->indexMap.insert( typename IndexMap::value_type( id, index ) );
	this->idVector[ index ] = id;

	if( index != 0 )
	{
	    moveUpPos( index ); 
	}


	sizeCheck();
	consistencyCheck();

	return id;
    }


    const bool isEmpty() const
    {
	return ( getSize() == 0 );
    }

    const Index getSize() const
    {
	return this->size;
    }


    void clear();


    void dump() const;

protected:

    void moveUp( const Index index )
    {
	const Index position( this->positionVector[index] );
	moveUpPos( position );
    }

    void moveDown( const Index index )
    {
	const Index position( this->positionVector[index] );
	moveDownPos( position );
    }

private:


    inline void moveUpPos( const Index position );
    inline void moveDownPos( const Index position );

    inline void sizeCheck() const;
    inline void consistencyCheck() const;

private:

    ItemVector    itemVector;
    IndexVector   heap;

    // maps itemVector index to heap position.
    IndexVector   positionVector;

    // map itemVector index to id.
    IDVector      idVector;

    // map id to itemVector index.
    IndexMap      indexMap;

    ID   idCounter;
    Index    size;

    PtrGreater< const Item* const > comp;

};




template < typename Item >
DynamicPriorityQueue< Item >::DynamicPriorityQueue()
    :
    size( 0 ),
    idCounter( 0 )
{
    ; // do nothing
}


template < typename Item >
void DynamicPriorityQueue< Item >::clear()
{
    this->itemVector.clear();
    this->heap.clear();
    this->positionVector.clear();
  
    this->size = 0;
}





template < typename Item >
void DynamicPriorityQueue< Item >::
movePos( Index pos )
{
    const Index index( this->heap[pos] );

    moveDownPos( pos );

    // If above moveDown() didn't move this item,
    // then we need to try moveUp() too.  If moveDown()
    // did work, nothing should be done.
    if( this->positionVector[index] == pos )
    {
	moveUpPos( pos );
    }
}


template < typename Item >
void DynamicPriorityQueue<Item>::moveUpPos( const Index position )
{
    const Index index( this->heap[position] );
    Item* const item( &getItemByIndex( index ) );

    Index pred( ( position - 1 ) / 2 );

    // first pass: do nothing if move up doesn't occur.
    Item* predItem( &getItemByIndex( this->heap[ pred ] ) );
    if( this->comp( item, predItem ) )
    {
	return;
    }

    // main loop
    Index pos( position );
    while( 1 )
    {
	this->heap[ pos ] = this->heap[ pred ];
	this->positionVector[ this->heap[ pred ] ] = pos;
	pos = pred;
      
	if( pred == 0 )
	{
	    break;
	}

	pred = ( pred - 1 ) / 2;

	predItem = &getItemByIndex( this->heap[ pred ] );

	if( this->comp( item, predItem ) )
	{
	    break;
	}
    }

    this->heap[pos] = index;
    this->positionVector[ index ] = pos;

    puts("u a");
    consistencyCheck();
    puts("u b");
}

// this is an optimized version.
template < typename Item >
void DynamicPriorityQueue< Item >::moveDownPos( const Index position )
{
    printf("dropping pos %d\n",position);
    assert( this->heap[ 0 ] <= heap.size() );
    const Index index( this->heap[position] );
    Item* const item( &getItemByIndex( index ) );

    Index succ( position * 2 + 1 );
 
    // first pass: simply return doing nothing if move down doesn't occur.
    if( succ >= heap.size() )
    {
	return;
    }

    if( this->comp( &getItemByIndex( this->heap[ succ ] ), 
		    &getItemByIndex( this->heap[ succ + 1 ] ) ) )
    {
	++succ;
    }
    
    if( this->comp( &getItemByIndex( this->heap[ succ ] ), item ) )
    {
	return;    // the going down does not occur.
    }
  
    Index pos( position );

    // main loop
    while( 1 )
    {
	// bring up the successor
	this->heap[ pos ] = this->heap[ succ ];
	this->positionVector[ this->heap[ succ ] ] = pos;
	pos = succ;

	// the next successor
	succ = succ * 2 + 1;

	if( succ >= heap.size() )
	{
	    break;
	}

	if( this->comp( &getItemByIndex( this->heap[ succ ] ), 
			&getItemByIndex( this->heap[ succ + 1 ] ) ) )
	{
	    ++succ;
	}
	
	// if the going down is finished, break.
	if( this->comp( &getItemByIndex( this->heap[ succ ] ), item ) )
	{
	    break;
	}
    }
    dump();
    puts("d a");
    printf("i %d p %d\n",index,pos);
    this->heap[pos] = index;
    this->positionVector[ index ] = pos;
    dump();

    consistencyCheck();
    puts("d b");
}


/* original version
   template < typename Item >
   void DynamicPriorityQueue< Item >::moveDown( Index index )
   {
   Index successor( index * 2 + 1 );

   if( successor < heap.size() - 1 && this->comp( this->heap[successor], this->heap[successor + 1] ) )
   {
   ++successor;
   }

   Item* item( this->heap[index] );
  
   while( successor < heap.size() && this->comp( item, this->heap[successor] ) )
   {
   this->heap[index] = this->heap[successor];
   this->positionVector[ this->heap[index] - theFirstItemPtr ] = index;
   index = successor;
   successor = index * 2 + 1;

   if( successor < heap.size() - 1 && 
   this->comp( this->heap[successor], this->heap[ successor + 1 ] ) )
   {
   ++successor;
   }
   }

   this->heap[index] = item;
   this->positionVector[ this->heap[index] - theFirstItemPtr ] = index;
   }
*/


template < typename Item >
void DynamicPriorityQueue< Item >::popItemByIndex( const Index index )
{
    puts("p a");

    dump();
    consistencyCheck();
    puts("p b");
    --this->size;
    
    // first, pop the item from the itemVector.
    this->itemVector[ index ] = this->itemVector[ getSize() ];
    this->itemVector.resize( getSize() );

    // update the idVector and the indexMap.
    printf("ind %d size %d\n",index,getSize());
    assert( index <= getSize() );
    const ID removedID( this->idVector[ index ] );
    const ID movedID( this->idVector[ getSize() ] );
    this->idVector[ index ] = movedID;
    this->idVector.resize( getSize() );

    for( Index i( 0 ); i < getSize(); ++i )
    {
	printf("%d ",idVector[i]);
    }

    printf("\nremoved %d moved %d\n", removedID, movedID );

    assert( this->indexMap.find( movedID ) != this->indexMap.end() );
    this->indexMap[ movedID ] = index;
    this->indexMap.erase( removedID );

    //
    // update the positionVector and the heap.
    //
    const Index removedPos( this->positionVector[ index ] );
    const Index movedPos( this->positionVector[ getSize() ] );
    printf("removedPos %d movedPos %d\n", removedPos, movedPos );
    
    // 1. swap positionVector[ end ] and positionVector[ index ]
    this->positionVector[ index ] = movedPos;
    this->heap[ movedPos ] = index;

    // 2. swap heap[ end ] and heap[ removed ].
    this->positionVector[ this->heap[ getSize() ] ] = removedPos;
    this->heap[ removedPos ] = this->heap[ getSize() ];

    // 3. discard the last.
    this->positionVector.resize( getSize() );
    this->heap.resize( getSize() );
    dump();
    if( getSize() != 0 )
    {
	movePos( removedPos );
    }

    assert( this->heap[ 0 ] <= getSize() );

    sizeCheck();
    consistencyCheck();
    puts("p c");
    dump();
}


template < typename Item >
void DynamicPriorityQueue< Item >::dump() const
{
    for( Index i( 0 ); i < getSize(); ++i )
    {
	printf("heap %d %d %d\n", i,heap[i],itemVector[heap[i]]);
    }
    for( Index i( 0 ); i < getSize(); ++i )
    {
	printf("pos %d %d\n", i,positionVector[i]);
    }
//    for( Index i( 0 ); i < getSize(); ++i )
//    {
//	printf("id %d\n",id);
//    }

}

template < typename Item >
void DynamicPriorityQueue< Item >::sizeCheck() const
{
#ifndef NDEBUG
    assert( this->itemVector.size() == getSize() );
    assert( this->heap.size() == getSize() );
    assert( this->positionVector.size() == getSize() );
    assert( this->idVector.size() == getSize() );
    assert( this->indexMap.size() == getSize() );
#endif /* NDEBUG */
}


template < typename Item >
void DynamicPriorityQueue< Item >::consistencyCheck() const
{
#ifndef NDEBUG
    // assert correct mapping between the heap and the positionVector.
    for( Index i( 0 ); i < getSize(); ++i )
    {
	assert( this->heap[ i ] <= getSize() );
	assert( this->positionVector[ i ] <= getSize() );
	assert( this->heap[ this->positionVector[i] ] == i );
    }


    // assert correct mapping between the indexMap and the idVector.
    for( Index i( 0 ); i < getSize(); ++i )
    {
	const ID id( this->idVector[i] );
	assert( id < this->idCounter );
	assert( this->indexMap.at( id ) == i );
    }

    // assert correct ordering of items in the heap.
    for( Index pos( 0 ); pos < getSize(); ++pos )
    {
	const Item* const item( &getItemByIndex( this->heap[ pos ] ) );

	const Index pred( ( pos - 1 ) / 2 );
	if( pred > 0 )
	{
	    const Item* const 
		predItem( &getItemByIndex( this->heap[ pred ] ) );

	    assert( this->comp( item, predItem ) );
	}

	const Index succ( pos * 2 + 1 );
	if( succ < getSize() )
	{
	    const Item* const 
		succItem( &getItemByIndex( this->heap[ succ ] ) );

	    assert( this->comp( succItem, item ) );
	}

    }

#endif /* NDEBUG */
}



#endif // __DYNAMICPRIORITYQUEUE_HPP



/*
  Do not modify
  $Author: shafi $
  $Revision: 2529 $
  $Date: 2005-11-19 01:36:40 -0800 (Sat, 19 Nov 2005) $
  $Locker$
*/





