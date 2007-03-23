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

    typedef typename ItemVector::size_type       Index;

    typedef std::vector< Index >  IndexVector;
    typedef std::vector< ID >     IDVector;

//    typedef std::tr1::unordered_map<const ID, Index> IndexMap;
    typedef std::map<const ID, Index> IndexMap;


    DynamicPriorityQueue();
  
    inline void move( const Index index );

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
	const Index index( this->size );
    
	++this->size;
    
//	    this->itemVector.resize( getSize() );
//	    this->heap.resize( getSize() );
//	    this->positionVector.resize( getSize() );
	
	this->itemVector.push_back( item );
	this->heap.push_back( index );
	this->positionVector.push_back( index );
	
	if( index != 0 )
	{
	    moveUpPos( index ); 
	}
	
	const ID id( this->idCounter );
	++this->idCounter;

	this->indexMap[ id ] = index;
	this->idVector.push_back( id );

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
move( Index index )
{
    //  assert( position < getSize() );
    const Index position( this->positionVector[index] );

    moveDownPos( position );

    // If above moveDown() didn't move this item,
    // then we need to try moveUp() too.  If moveDown()
    // did work, nothing should be done.
    if( this->positionVector[index] == position )
    {
	moveUpPos( position );
    }
}


template < typename Item >
void DynamicPriorityQueue<Item>::moveUpPos( const Index position )
{
    const Index index( this->heap[position] );
    Item* const item( &getItemByIndex( index ) );

    Index predecessor( ( position - 1 ) / 2 );

    // first pass: do nothing if move up doesn't occur.
    Item* predItem( &getItemByIndex( this->heap[predecessor] ) );
    if( this->comp( item, predItem ) )
    {
	return;
    }

    // main loop
    Index pos( position );
    while( 1 )
    {
	this->heap[pos] = predecessor;
	this->positionVector[ this->heap[predecessor] ] = pos;
	pos = predecessor;
      
	if( predecessor == 0 )
	{
	    break;
	}

	predecessor = ( predecessor - 1 ) / 2;

	predItem = &getItemByIndex( this->heap[predecessor] );

	if( this->comp( item, predItem ) )
	{
	    break;
	}
    }

    this->heap[pos] = index;
    this->positionVector[ index ] = position;
}

// this is an optimized version.
template < typename Item >
void DynamicPriorityQueue< Item >::moveDownPos( Index position )
{
    Index index( this->heap[position] );
    Item* const item( &getItemByIndex( index ) );
    Index successor( position * 2 + 1);
 

    // first pass: simply return doing nothing if move down doesn't occur.
    if( successor < getSize() - 1 )
    {
	if( this->comp( &getItemByIndex( this->heap[ successor ] ), 
			&getItemByIndex( this->heap[ successor + 1 ] ) ) )
	{
	    ++successor;
	}
    }
    else if( successor >= getSize() )
    {
	return;
    }
  
    Item* succItem( &getItemByIndex( this->heap[ successor ] ) );
    if( this->comp( succItem, item ) )
    {
	return;    // if the going down does not occur, return doing nothing.
    }

    // main loop
    while( 1 )
    {
	// bring up the successor
	this->heap[position] = successor;
	this->positionVector[ this->heap[ successor ] ] = position;
	position = successor;

	// the next successor
	successor = successor * 2 + 1;

	if( successor < getSize() - 1 )
	{
	    if( this->comp( &getItemByIndex( this->heap[ successor ] ), 
			    &getItemByIndex( this->heap[ successor + 1 ] ) ) )
	    {
		++successor;
	    }
	}
	else if( successor >= getSize() )
	{
	    break;
	}

	succItem = &getItemByIndex( this->heap[ successor ] );

	// if the going down is finished, break.
	if( this->comp( succItem, item ) )
	{
	    break;
	}
    }

    this->heap[position] = index;
    this->positionVector[ index ] = position;
}


/* original version
   template < typename Item >
   void DynamicPriorityQueue< Item >::moveDown( Index index )
   {
   Index successor( index * 2 + 1 );

   if( successor < getSize() - 1 && this->comp( this->heap[successor], this->heap[successor + 1] ) )
   {
   ++successor;
   }

   Item* item( this->heap[index] );
  
   while( successor < getSize() && this->comp( item, this->heap[successor] ) )
   {
   this->heap[index] = this->heap[successor];
   this->positionVector[ this->heap[index] - theFirstItemPtr ] = index;
   index = successor;
   successor = index * 2 + 1;

   if( successor < getSize() - 1 && 
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
    --this->size;


    // first, pop the item from the itemVector.
    this->itemVector[ index ] = this->itemVector[ getSize() ];
    this->itemVector.resize( getSize() );

    // update the idVector.
    const ID movedID( this->idVector[ getSize() ] );
    const ID removedID( this->idVector[ index ] );
    this->idVector[ index ] = this->idVector[ getSize() ];
    this->idVector.resize( getSize() );

    // update the positionVector
    const Index pos( this->positionVector[ index ] );
    const Index movedPos( this->positionVector[ getSize() ] );
    this->positionVector[ index ] = this->positionVector[ getSize() ];
    this->positionVector.resize( getSize() );


    // 
    this->heap[pos] = movedPos;
    moveDownPos( pos );

    this->heap.resize( getSize() );



    printf("i p %d %d\n",index,pos);


	
    this->indexMap.erase( removedID );
    this->indexMap[ movedID ] = index;


    if( getSize() != 0 )
    {

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





