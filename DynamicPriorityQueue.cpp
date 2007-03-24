
#include "DynamicPriorityQueue.hpp"


#ifdef DPQ_TEST

#include <iostream>

int main()
{
  DynamicPriorityQueue<int> dpq;


  for( int i( 4 ); i != 0  ; i-=2 )
    {
	dpq.pushItem( i );
	std::cout << i << ' ' << dpq.getTopItem() << std::endl;
    }

  for( int i( 5 ); i != -1  ; i-=2 )
    {
	dpq.pushItem( i );
	std::cout << i << ' ' << dpq.getTopItem() << std::endl;
    }

//  dpq.dump();


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
//	dpq.dump();
    }


}


#endif /* DPQ_TEST */
