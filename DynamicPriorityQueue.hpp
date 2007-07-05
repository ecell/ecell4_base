//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//        This file is part of E-Cell Simulation Environment package
//
//                Copyright (C) 1996-2007 Keio University
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
// modified by Mozoyoshi Koizumi
//

#ifndef __DYNAMICPRIORITYQUEUE_HPP
#define __DYNAMICPRIORITYQUEUE_HPP

#include "config.h"

#include <functional>
#include <vector>
#include <algorithm>
#include <stdexcept>

//#define HAVE_TR1_UNORDERED_MAP

#if HAVE_UNORDERED_MAP
#include <unordered_map>
#elif HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
#else
#include <map>
#endif /* HAVE_UNORDERED_MAP */


//namespace libecs
//{


class PersistentIDPolicy
{
    
public:
    
    typedef long long unsigned int ID;
    
    typedef std::vector< ID >      IDVector;
    typedef IDVector::size_type    Index;

#if HAVE_UNORDERED_MAP || HAVE_TR1_UNORDERED_MAP

    class IDHasher
        : 
        public std::unary_function<ID, std::size_t>
    {

    public:

        std::size_t operator()( ID value ) const
        {
            return static_cast<std::size_t>( value ) ^
                static_cast<std::size_t>( value >> ( sizeof( ID ) * 8 / 2 ) );
        }

    };

#endif // HAVE_UNORDERED_MAP || HAVE_TR1_UNORDERED_MAP

#if HAVE_UNORDERED_MAP
    typedef std::unordered_map<const ID, Index, IDHasher> IndexMap;
#elif HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<const ID, Index, IDHasher> IndexMap;
#else 
    typedef std::map<const ID, Index> IndexMap;
#endif


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
        IndexMap::const_iterator i = this->indexMap.find( id );

        if( i == this->indexMap.end() )
        {
            throw std::out_of_range( "PersistentIDPolicy::getIndex():"
				     " Key not found." );
        }

        return (*i).second;
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

    const bool check( const Index size ) const
    {
        if( this->idVector.size() != size )
        {
            return false;
        }

        if( this->indexMap.size() != size )
        {
            return false;
        }

        // assert correct mapping between the indexMap and the idVector.
        for( Index i( 0 ); i < size; ++i )
        {
            const ID id( this->idVector[i] );

            if (id >= this->idCounter)
            {
                return false;
            }

            IndexMap::const_iterator iter = this->indexMap.find( id );
            if (iter == this->indexMap.end())
            {
                return false;
            }

            if ((*iter).second != i)
            {
                return false;
            }
        }

        return true;
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

    const bool check( const Index size ) const
    {
        return true;
    }


};


/**
   Dynamic priority queue for items of type Item.

   When IDPolicy template parameter is PersistentIDPolicy, IDs assigned
   to pushed items are persistent for the life time of this priority
   queue.

   When VolatileIDPolicy template parameter is used as the IDPolicy,
   IDs are valid only until the next call of pop or push methods.
   However, VolatileIDPolicy saves some memory and eliminates the
   overhead incurred in pop/push methods.
*/

template < typename Item, class IDPolicy = PersistentIDPolicy >
class DynamicPriorityQueue
    :
    private IDPolicy
{

public:

    typedef typename IDPolicy::ID    ID;
    typedef typename IDPolicy::Index Index;

    typedef std::vector< Item >      ItemVector;
    typedef std::vector< Index >     IndexVector;


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

    Item& getTop()
    {
        return this->itemVector[ getTopIndex() ];
    }

    const Item& getTop() const
    {
        return this->itemVector[ getTopIndex() ];
    }

    Item& get( const ID id )
    {
        return this->itemVector[ getIndex( id ) ];
    }

    const Item& get( const ID id ) const
    {
        return this->itemVector[ getIndex( id ) ];
    }

    ID getTopID() const
    {
        return IDPolicy::getIDByIndex( getTopIndex() );
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

    const Item& getByIndex( const Index index ) const
    {
        return this->itemVector[ index ];
    }

    const Index getTopIndex() const 
    {
        return this->heap[0];
    }

    Item& operator[]( const ID id )
    {
        return get( id );
    }

    const Item& operator[]( const ID id ) const
    {
        return get( id );
    }


    // self-diagnostic methods
    const bool check() const; // check all
    const bool checkSize() const;
    const bool checkPositionMapping() const;
    const bool checkHeap() const;
    const bool checkIDPolicy() const;

    void dump() const;

private:

    inline void popByIndex( const Index index );

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


    inline void moveUpPos( const Index position, const Index start = 0 );
    inline void moveDownPos( const Index position );

    inline void moveUpPosNoCheck( const Index position, 
                                       const Index start = 0 );
    inline void moveDownPosNoCheck( const Index position ); 

private:

    ItemVector    itemVector;
    IndexVector   heap;

    // maps itemVector index to heap position.
    IndexVector   positionVector;


    std::less_equal< const Item > comp;
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
        if( ! this->comp( item, this->itemVector[ this->heap[ succ ] ] ) ||
            ( succ + 1 < size && 
              ! this->comp( item, 
                            this->itemVector[ this->heap[ succ + 1 ] ] ) ) )
        {
            moveDownPosNoCheck( pos );
            return;
        }
    }

    if( pos <= 0 ) // already on top
    {
        return;
    }

    moveUpPos( pos );
}

template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::moveUpPos( const Index position, 
                                                        const Index start )
{

    const Index index( this->heap[ position ] );
    const Item& item( this->itemVector[ index ] );

    const Index pred( ( position - 1 ) / 2 );
    const Index predIndex( this->heap[ pred ] );

    if( ! this->comp( this->itemVector[ predIndex ], item ) )
    {
        this->moveUpPosNoCheck( position, start );
    }
}


template < typename Item, class IDPolicy >
void 
DynamicPriorityQueue< Item, IDPolicy >::moveDownPos( const Index position )
{
    const Index index( this->heap[ position ] );
    const Item& item( this->itemVector[ index ] );

    const Index size( getSize() );

    const Index succ( 2 * position + 1 );
    if( succ < size )
    {
        if( ! this->comp( item, this->itemVector[ this->heap[ succ ] ] ) ||
            ( succ + 1 < size && 
              ! this->comp( item, 
                            this->itemVector[ this->heap[ succ + 1 ] ] ) ) )
        {
            this->moveDownPosNoCheck( position );
        }
    }
}

template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >:: 
moveUpPosNoCheck( const Index position,
                  const Index start )
{
    const Index index( this->heap[ position ] );
    const Item& item( this->itemVector[ index ] );

    if( position <= start )
    {
        return;
    }

    Index pos( position );
    Index pred( ( pos - 1 ) / 2 );
    Index predIndex( this->heap[ pred ] );

    do
    {
        this->heap[ pos ] = predIndex;
        this->positionVector[ predIndex ] = pos;
        pos = pred;

        if( pos <= start )
        {
            break;
        }

        pred = ( pos - 1 ) / 2;
        predIndex = this->heap[ pred ];

    } while( ! this->comp( this->itemVector[ predIndex ], item ) );


    this->heap[ pos ] = index;
    this->positionVector[ index ] = pos;
}


template < typename Item, class IDPolicy >
void 
DynamicPriorityQueue< Item, IDPolicy >::
moveDownPosNoCheck( const Index position )
{
    const Index index( this->heap[ position ] );

    const Index size( this->getSize() );
    
    Index succ( 2 * position + 1 );
    Index pos( position );
    while( succ < size )
    {
        const Index rightPos( succ + 1 );
        if( rightPos < size && 
            ! this->comp( this->itemVector[ this->heap[ succ ] ],
                          this->itemVector[ this->heap[ rightPos ] ] ) )
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

    // assert( check() );

    return id;
}


template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::popByIndex( const Index index )
{
    // 1. pop the item from the itemVector.
    this->itemVector[ index ] = this->itemVector.back();
    this->itemVector.pop_back();

    // 2. update index<->ID mapping.
    IDPolicy::pop( index );

    const Index removedPos( this->positionVector[ index ] );
    const Index movedPos( this->positionVector.back() );

    // 3. swap positionVector[ end ] and positionVector[ index ]
    this->positionVector[ index ] = movedPos;
    this->heap[ movedPos ] = index;

    // 4. if heap[ end ] and heap[ removed ] do not overlap,
    //    swap these, pop back, and update the heap.
    if( removedPos != this->heap.size() - 1 )
    {
        this->heap[ removedPos ] = this->heap.back();
        this->positionVector[ this->heap.back() ] = removedPos;

        this->positionVector.pop_back();
        this->heap.pop_back();

        movePos( removedPos );
    }
    else  // if heap[ end ] and heap[ removed ] are the same, simply pop back.
    {
        this->positionVector.pop_back();
        this->heap.pop_back();
    }

    //    assert( check() );
}



template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::replaceTop( const Item& item )
{
    this->itemVector[ this->heap[0] ] = item;
    moveTop();
    
    // assert( check() );
}

template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::
replace( const ID id, const Item& item )
{
    const Index index( getIndex( id ) );
    this->itemVector[ index ] = item;
    move( index );
    
    // assert( check() );
}



#include <iostream>

template < typename Item, class IDPolicy >
void DynamicPriorityQueue< Item, IDPolicy >::dump() const
{
    for( Index i( 0 ); i < heap.size(); ++i )
    {
        std::cerr << "heap\t" << i << " " << heap[i] << std::endl;
    }
    for( Index i( 0 ); i < positionVector.size(); ++i )
    {
        std::cerr << "pos\t" << i << " " << positionVector[i] << std::endl;
    }
}


template < typename Item, class IDPolicy >
const bool DynamicPriorityQueue< Item, IDPolicy >::check() const
{
    bool result( true );

    result = result && this->checkSize();
    result = result && this->checkPositionMapping();
    result = result && this->checkHeap();
    result = result && this->checkIDPolicy();

    return result;
}


template < typename Item, class IDPolicy >
const bool DynamicPriorityQueue< Item, IDPolicy >::checkSize() const
{
    bool result( true );

    // check sizes of data structures.
    result = result && this->itemVector.size() == getSize();
    result = result && this->heap.size() == getSize();
    result = result && this->positionVector.size() == getSize();

    return result;
}


template < typename Item, class IDPolicy >
const bool DynamicPriorityQueue< Item, IDPolicy >::checkPositionMapping() const
{
    bool result( true );

    // assert correct mapping between the heap and the positionVector.
    for( Index i( 0 ); i < getSize(); ++i )
    {
        result = result && this->heap[ i ] < getSize();
        result = result && this->positionVector[ i ] < getSize();
        result = result && this->heap[ this->positionVector[i] ] == i;
    }

    return result;
}

template < typename Item, class IDPolicy >
const bool DynamicPriorityQueue< Item, IDPolicy >::checkHeap() const
{
    bool result( true );

    // assert correct ordering of items in the heap.

    for( Index pos( 0 ); pos < getSize(); ++pos )
    {
        const Item& item( this->itemVector[ this->heap[ pos ] ] );

        const Index succ( pos * 2 + 1 );
        if( succ < getSize() )
        {
            result = result && 
                this->comp( item, this->itemVector[ this->heap[ succ ] ] );

            const Index rightPos( succ + 1 );
            if( rightPos < getSize() )
            {
                result = result &&  
                    this->comp( item, 
                                this->itemVector[ this->heap[ rightPos ] ] );
            }
        }

    }

    return result;
}


template < typename Item, class IDPolicy >
const bool DynamicPriorityQueue< Item, IDPolicy >::checkIDPolicy() const
{
    bool result( true );

    result = result && IDPolicy::check( getSize() );
    return result;
}



//} // namespace libecs

#endif // __DYNAMICPRIORITYQUEUE_HPP
