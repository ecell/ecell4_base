
#include "DynamicPriorityQueue.hpp"


#ifdef DPQ_TEST

#include <iostream>

int main()
{
  DynamicPriorityQueue<int> dpq;


  for( int i( 10 ); i != 0  ; --i )
    {
	dpq.pushItem( i );
	std::cout << i << ' ' << dpq.getTopItem() << std::endl;
    }

//  dpq.popItem( 3 );
//  dpq.popItem( 5 );
//  std::cerr << dpq.getTopItem() << std::endl;
//  dpq.popTop();
//  std::cerr << dpq.getTopItem() << std::endl;

//  DynamicPriorityQueue<int> copy( dpq );
  while( ! dpq.isEmpty() )
    {
	std::cerr << dpq.getSize() << ' ' << dpq.getTopItem() << std::endl;
	dpq.popTop();
    }


}


#endif /* DPQ_TEST */
