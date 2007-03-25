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

    typedef typename ItemVector::difference_type       Index;

    typedef std::vector< Index >   IndexVector;
    typedef std::vector< ID >      IDVector;

//    typedef std::tr1::unordered_map<const ID, Index> IndexMap;
    typedef std::map<const ID, Index> IndexMap;


    DynamicPriorityQueue();
  
    void move( const Index index )
    {
	const Index pos( this->positionVector[ index ] );
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
	return this->heap[0];
    }

    const Item& getTopItem() const
    {
	return this->itemVector[ getTopIndex() ];
    }

    Item& getTopItem()
    {
	return this->itemVector[ getTopIndex() ];
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
	return this->getItemByIndex( this->indexMap.at( id ) );
    }

    Item& getItem( const ID id )
    {
	return this->getItemByIndex( this->indexMap.at( id ) );
    }

    void popTop()
    {
	return popItemByIndex( getTopIndex() );
    }

    void popItem( const ID id )
    {
	return popItemByIndex( getIndex( id ) );
    }

    inline void popItemByIndex( const Index index );

    inline const ID pushItem( const Item& item );

    const bool isEmpty() const
    {
	return this->itemVector.empty();
    }

    const Index getSize() const
    {
	return this->itemVector.size();
    }

    void clear();

    void dump() const;

protected:

    void moveUp( const Index index )
    {
	const Index position( this->positionVector[ index ] );
	moveUpPos( position );
    }

    void moveDown( const Index index )
    {
	const Index position( this->positionVector[ index ] );
	moveDownPos( position );
    }

private:

    inline void moveUpPos( const Index position, const Index start = 0 );
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

    PtrGreater< const Item* const > comp;
};




template < typename Item >
DynamicPriorityQueue< Item >::DynamicPriorityQueue()
    :
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
    this->idVector.clear();
    this->indexMap.clear();
}


template < typename Item >
void DynamicPriorityQueue< Item >::
movePos( const Index pos )
{
    const Index index( this->heap[ pos ] );
    const Item* const item( &getItemByIndex( index ) );

    const Index size( getSize() );

    const Index succ( 2 * pos + 1 );
    if( succ < size )
    {
	if( this->comp( item, &getItemByIndex( this->heap[ succ ] ) ) ||
	    ( succ + 1 < size && 
	      this->comp( item, &getItemByIndex( this->heap[ succ + 1 ] ) ) ) )
	{
	    moveDownPos( pos );
	    return;
	}
    }

    const Index pred( ( pos - 1 ) / 2 );
    if( pred >= 0  && 
	this->comp( &getItemByIndex( this->heap[ pred ] ), item ) )
    {
	moveUpPos( pos );
    }
}

template < typename Item >
void DynamicPriorityQueue<Item>::moveUpPos( const Index position, 
					    const Index start )
{
    const Index index( this->heap[ position ] );
    const Item* const item( &getItemByIndex( index ) );

    Index pos( position );
    while( pos > start )
    {
	const Index pred( ( pos - 1 ) / 2 );
	const Index predIndex( this->heap[ pred ] );
	if( this->comp( item, &getItemByIndex( predIndex ) ) )
	{
	    break;
	}

	this->heap[ pos ] = predIndex;
	this->positionVector[ predIndex ] = pos;
	pos = pred;
    }

    this->heap[ pos ] = index;
    this->positionVector[ index ] = pos;
}


template < typename Item >
void DynamicPriorityQueue< Item >::moveDownPos( const Index position )
{
    const Index index( this->heap[ position ] );
    const Item* const item( &getItemByIndex( index ) );

    const Index size( getSize() );
    
    Index succ( 2 * position + 1 );
    Index pos( position );
    while( succ < size )
    {
	const Index rightPos( succ + 1 );
	if( rightPos < size && 
	    this->comp( &getItemByIndex( this->heap[ succ ] ),
			&getItemByIndex( this->heap[ rightPos ] ) ) )
	{
	    succ = rightPos;
	}

	this->heap[ pos ] = this->heap[ succ ];
	this->positionVector[ this->heap[ pos ] ] = pos;
	pos = succ;
	succ = 2 * pos + 1;
    }

    this->heap[ pos ] = index;
    this->positionVector[ index ] = pos;

    moveUpPos( pos, position );
}


template < typename Item >
const typename DynamicPriorityQueue<Item>::ID
DynamicPriorityQueue<Item>::pushItem( const Item& item )
{
    const Index index( getSize() );
    
    this->itemVector.push_back( item );
    // index == pos at this time.
    this->heap.push_back( index );
    this->positionVector.push_back( index );

    const ID id( this->idCounter );
    ++this->idCounter;

    this->indexMap.insert( typename IndexMap::value_type( id, index ) );
    this->idVector.push_back( id );

    moveUpPos( index ); 

//    sizeCheck();
//    consistencyCheck();

    return id;
}


template < typename Item >
void DynamicPriorityQueue< Item >::popItemByIndex( const Index index )
{
    // first, pop the item from the itemVector.
    this->itemVector[ index ] = this->itemVector.back();
    this->itemVector.pop_back();

    // update the idVector and the indexMap.
    const ID removedID( this->idVector[ index ] );
    const ID movedID( this->idVector.back() );
    this->idVector[ index ] = movedID;
    this->idVector.pop_back();

    this->indexMap[ movedID ] = index;
    this->indexMap.erase( removedID );

    //
    // update the positionVector and the heap.
    //
    const Index removedPos( this->positionVector[ index ] );
    const Index movedPos( this->positionVector.back() );
    
    // 1. swap positionVector[ end ] and positionVector[ index ]
    this->positionVector[ index ] = movedPos;
    this->heap[ movedPos ] = index;

    // 2. swap heap[ end ] and heap[ removed ].
    this->positionVector[ this->heap.back() ] = removedPos;
    this->heap[ removedPos ] = this->heap.back();

    // 3. discard the last.
    this->positionVector.pop_back();
    this->heap.pop_back();

    movePos( removedPos );

//    sizeCheck();
//    consistencyCheck();
}


template < typename Item >
void DynamicPriorityQueue< Item >::dump() const
{
    for( Index i( 0 ); i < heap.size(); ++i )
    {
	printf("heap %d %d %d\n", i,heap[i],itemVector[heap[i]]);
    }
    for( Index i( 0 ); i < positionVector.size(); ++i )
    {
	printf("pos %d %d\n", i,positionVector[i]);
    }
}

template < typename Item >
void DynamicPriorityQueue< Item >::sizeCheck() const
{
    assert( this->itemVector.size() == getSize());
    assert( this->heap.size() == getSize());
    assert( this->positionVector.size() == getSize());
    assert( this->idVector.size() == getSize());
    assert( this->indexMap.size() == getSize());
}


template < typename Item >
void DynamicPriorityQueue< Item >::consistencyCheck() const
{
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

	const Index succ( pos * 2 + 1 );
	if( succ < getSize() )
	{
	    assert( this->comp( &getItemByIndex( this->heap[succ] ), item ) );

	    if( succ+1 < getSize() )
	    {
		assert( this->comp( &getItemByIndex( this->heap[succ+1] ),
				    item ));
	    }
	}

    }
}



#endif // __DYNAMICPRIORITYQUEUE_HPP



/*
  Do not modify
  $Author: shafi $
  $Revision: 2529 $
  $Date: 2005-11-19 01:36:40 -0800 (Sat, 19 Nov 2005) $
  $Locker$
*/





