
#include "DynamicPriorityQueue.hpp"


#ifdef DPQ_TEST

#include <iostream>

int main()
{
  DynamicPriorityQueue<int> dpq;


  for( int i( 100000 ); i != 0  ; i-=2 )
    {
	dpq.pushItem( i );
//	std::cout << i << ' ' << dpq.getTopItem() << std::endl;
    }

  for( int i( 100001 ); i != -1  ; i-=2 )
    {
	int id = dpq.pushItem( i );
//	std::cout << i << ' ' << dpq.getTopItem() << ' ' << id << std::endl;
    }

//  dpq.dump();

  std::cout << "popping " << dpq.getItem( 0 ) << std::endl;
  dpq.popItem( 0 );

  std::cout << "popping " << dpq.getItem( 9 ) << std::endl;
  dpq.popItem( 9 );
//  dpq.popItem( 5 );
//  std::cerr << dpq.getTopItem() << std::endl;
//  dpq.popTop();
//  std::cerr << dpq.getTopItem() << std::endl;

//  DynamicPriorityQueue<int> copy( dpq );
  while( ! dpq.isEmpty() )
    {
//	std::cout << "pop: " << dpq.getSize() << ' ' << dpq.getTopItem() << std::endl;
	dpq.popTop();
//	dpq.dump();
    }


}


#endif /* DPQ_TEST */
