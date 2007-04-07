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
// written by Koichi Takahashi based on the initial version by Eiichiro Adachi.
//

#ifndef __DYNAMICPRIORITYQUEUE_HPP
#define __DYNAMICPRIORITYQUEUE_HPP
//#include <assert.h>
#include <functional>
#include <vector>
#include <algorithm>
#include <map>
//#include <tr1/unordered_map>



class PersistentIDPolicy
{

public:

    typedef long long unsigned int ID;

    typedef std::vector< ID >      IDVector;
    typedef IDVector::size_type    Index;

//    typedef std::tr1::unordered_map<const ID, Index> IndexMap;
    typedef std::map<const ID, Index> IndexMap;

    PersistentIDPolicy()
	:
	idCounter( 0 )
    {
	; // do nothing
    }
    
    void reset()
    {
	idCounter = 0;
    }

    void clear()
    {
	this->idVector.clear();
	this->indexMap.clear();
    }

    const Index getIndex( const ID id ) const
    {
	return this->indexMap.at( id );
    }

    const ID getIDByIndex( const Index index ) const
    {
	return this->idVector[ index ];
    }

    const ID push( const Index index )
    {
	const ID id( this->idCounter );
	++this->idCounter;

	this->indexMap.insert( IndexMap::value_type( id, index ) );
	this->idVector.push_back( id );

	return id;
    }

    void pop( const Index index )
    {
	// update the idVector and the indexMap.
	const ID removedID( this->idVector[ index ] );
	const ID movedID( this->idVector.back() );
	this->idVector[ index ] = movedID;
	this->idVector.pop_back();
	
	this->indexMap[ movedID ] = index;
	this->indexMap.erase( removedID );
    }

    const bool checkConsistency( const Index size ) const
    {
	bool result( true );

	result &= this->idVector.size() == size;
	result &= this->indexMap.size() == size;

	// assert correct mapping between the indexMap and the idVector.
	for( Index i( 0 ); i < size; ++i )
	{
	    const ID id( this->idVector[i] );
	    result &= id < this->idCounter;
	    result &= this->indexMap.at( id ) == i;
	}

	return result;
    }


private:

    // map itemVector index to id.
    IDVector      idVector;
    // map id to itemVector index.
    IndexMap      indexMap;

    ID   idCounter;

};

class VolatileIDPolicy
{
public:


    typedef size_t    Index;
    typedef Index     ID;

    void reset()
    {
	; // do nothing
    }

    void clear()
    {
	; // do nothing
    }

    const Index getIndex( const ID id ) const
    {
	return id;
    }

    const ID getIDByIndex( const Index index ) const
    {
	return index;
    }

    const ID push( const Index index )
    {
	return index;
    }

    void pop( const Index index )
    {
	; // do nothing
    }

    const bool checkConsistency( const Index size ) const
    {
	return true;
    }


};

/**
   Dynamic priority queue for items of type Item.

   When IDPolicy is PersistentIDPolicy, IDs given to
   pushed items are persistent for the life time of this
   priority queue.

   When VolatileIDPolicy is used as the IDPolicy, IDs
   are valid only until the next call or pop or push methods.
   However, VolatileIDPolicy saves some memory and eliminates
   the overhead incurred in pop/push methods.

*/

template < typename Item, class IDPolicy = PersistentIDPolicy >
class DynamicPriorityQueue
    :
    private IDPolicy
{
  

public:

    typedef std::vector< Item >    ItemVector;

    typedef typename IDPolicy::ID ID;
    typedef typename IDPolicy::Index Index;

    typedef std::vector< Index >   IndexVector;


    DynamicPriorityQueue();
  
    const bool isEmpty() const
    {
	return this->itemVector.empty();
    }

    const Index getSize() const
    {
	return this->itemVector.size();
    }

    void clear();

    const Item& getTop() const
    {
	return this->itemVector[ getTopIndex() ];
    }

    Item& getTop()
    {
	return this->itemVector[ getTopIndex() ];
    }

    ID getTopID() const
    {
	return IDPolicy::getIDByIndex( getTopIndex() );
    }

    const Item& operator[]( const ID id ) const
    {
	return this->itemVector[ getIndex( id ) ];
    }

    Item& operator[]( const ID id )
    {
	return this->itemVector[ getIndex( id ) ];
    }

    void popTop()
    {
	popByIndex( getTopIndex() );
    }

    void pop( const ID id )
    {
	popByIndex( getIndex( id ) );
    }

    void replaceTop( const Item& item );

    void replace( const ID id, const Item& item );

    inline const ID push( const Item& item );

    void dump() const;


    // self-diagnostic method
    const bool checkConsistency() const;

protected:

    inline void popByIndex( const Index index );

    Item& getByIndex( const Index index )
    {
	return this->itemVector[ index ];
    }

    const Index getTopIndex() const 
    {
	return this->heap[0];
    }

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

private:

    ItemVector    itemVector;
    IndexVector   heap;

    // maps itemVector index to heap position.
    IndexVector   positionVector;


    static std::less_equal< const Item > comp;
};




template < typename Item, class IDPolicy >
DynamicPriorityQueue< Item, IDPolicy >::DynamicPriorityQueue()
{
    ; // do nothing
}


template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::clear()
{
    this->itemVector.clear();
    this->heap.clear();
    this->positionVector.clear();
    IDPolicy::clear();
}


template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::
movePos( const Index pos )
{
    const Index index( this->heap[ pos ] );
    const Item& item( this->itemVector[ index ] );

    const Index size( getSize() );

    const Index succ( 2 * pos + 1 );
    if( succ < size )
    {
	if( this->comp( this->itemVector[ this->heap[ succ ] ], item ) ||
	    ( succ + 1 < size && 
	      this->comp( this->itemVector[ this->heap[ succ + 1 ] ], 
			  item ) ) )
	{
	    moveDownPos( pos );
	    return;
	}
    }

    if( pos <= 0 )
    {
	return;
    }

    const Index pred( ( pos - 1 ) / 2 );
    if( pred >= 0  && 
	this->comp( item, this->itemVector[ this->heap[ pred ] ] ) )
    {
	moveUpPos( pos );
	return;
    }
}

template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::moveUpPos( const Index position, 
							const Index start )
{
    const Index index( this->heap[ position ] );
    const Item& item( this->itemVector[ index ] );

    Index pos( position );
    while( pos > start )
    {
	const Index pred( ( pos - 1 ) / 2 );
	const Index predIndex( this->heap[ pred ] );
	if( this->comp( this->itemVector[ predIndex ], item ) )
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


template < typename Item, class IDPolicy >
void 
DynamicPriorityQueue< Item, IDPolicy >::moveDownPos( const Index position )
{
    const Index index( this->heap[ position ] );
    const Item& item( this->itemVector[ index ] );

    const Index size( getSize() );
    
    Index succ( 2 * position + 1 );
    Index pos( position );
    while( succ < size )
    {
	const Index rightPos( succ + 1 );
	if( rightPos < size && 
	    this->comp( this->itemVector[ this->heap[ rightPos ] ],
			this->itemVector[ this->heap[ succ ] ] ) )
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


template < typename Item, class IDPolicy >
const typename DynamicPriorityQueue< Item, IDPolicy >::ID
DynamicPriorityQueue< Item, IDPolicy >::push( const Item& item )
{
    const Index index( getSize() );
    
    this->itemVector.push_back( item );
    // index == pos at this time.
    this->heap.push_back( index );
    this->positionVector.push_back( index );

    const ID id( IDPolicy::push( index ) );

    moveUpPos( index ); 

//    assert( checkConsistency() );

    return id;
}


template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::popByIndex( const Index index )
{
    // first, pop the item from the itemVector.
    this->itemVector[ index ] = this->itemVector.back();
    this->itemVector.pop_back();


    // update index<->ID mapping.
    IDPolicy::pop( index );

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

//    assert( checkConsistency() );
}



template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::replaceTop( const Item& item )
{
    this->itemVector[ this->heap[0] ] = item;
    moveTop();
    
//    assert( checkConsistency() );
}

template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::
replace( const ID id, const Item& item )
{
    const Index index( getIndex( id ) );
    this->itemVector[ index ] = item;
    move( index );
    
//    assert( checkConsistency() );
}


template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::dump() const
{
    for( Index i( 0 ); i < heap.size(); ++i )
    {
	printf( "heap %d %d %d\n", 
		i, heap[i], this->itemVector[ this->heap[i] ] );
    }
    for( Index i( 0 ); i < positionVector.size(); ++i )
    {
	printf( "pos %d %d\n", 
		i, positionVector[i] );
    }
}


template < typename Item, class IDPolicy >
const bool DynamicPriorityQueue< Item, IDPolicy >::checkConsistency() const
{
    bool result( true );

    // check sizes of data structures.
    result &= this->itemVector.size() == getSize();
    result &= this->heap.size() == getSize();
    result &= this->positionVector.size() == getSize();

    // assert correct mapping between the heap and the positionVector.
    for( Index i( 0 ); i < getSize(); ++i )
    {
	result &= this->heap[ i ] <= getSize();
	result &= this->positionVector[ i ] <= getSize();
	result &= this->heap[ this->positionVector[i] ] == i;
    }

    // assert correct ordering of items in the heap.

    for( Index pos( 0 ); pos < getSize(); ++pos )
    {
	const Item& item( this->itemVector[ this->heap[ pos ] ] );

	const Index succ( pos * 2 + 1 );
	if( succ < getSize() )
	{
	    result &= this->comp( item, 
				  this->itemVector[ this->heap[ succ ] ] );

	    const Index rightPos( succ + 1 );
	    if( rightPos < getSize() )
	    {
		result &= 
		    this->comp( item, 
				this->itemVector[ this->heap[ rightPos ] ] );
	    }
	}

    }

    result &= IDPolicy::checkConsistency( getSize() );

    return result;
}



#endif // __DYNAMICPRIORITYQUEUE_HPP
