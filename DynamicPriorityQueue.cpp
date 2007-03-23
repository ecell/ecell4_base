
#include "DynamicPriorityQueue.hpp"


#ifdef DPQ_TEST

#include <iostream>

int main()
{
  DynamicPriorityQueue<int> dpq;


  for( int i( 100 ); i != 0  ; --i )
    //  for( int i( 0 ); i < 10000  ; ++i )
    {
	std::cout << i << std::endl;
	dpq.pushItem( i );
    }

  std::cerr << dpq.getTopItem() << std::endl;
  dpq.popTop();
  std::cerr << dpq.getTopItem() << std::endl;

/*
  for( int i( 0 ); i < 10000  ; ++i )
    {
      dpq.changeOneKey( dpq.topIndex(), i+1000 );
    }
*/

//  DynamicPriorityQueue<int> copy( dpq );
  while( ! dpq.isEmpty() )
    {
	std::cerr << dpq.getSize() << ' ' << dpq.getTopItem() << std::endl;
	dpq.getTopItem() = 10000;
	dpq.moveDown( dpq.getTopIndex() );
	dpq.popTop();
    }


}


#endif /* DPQ_TEST */
