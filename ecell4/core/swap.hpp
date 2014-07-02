#ifndef __ECELL4_UTILS_SWAP_HPP
#define __ECELL4_UTILS_SWAP_HPP

#include <cstring>

namespace ecell4
{

template<typename T>
void blit_swap(T& x, T& y)
{
    if (&x == &y)
        return;
    struct blob { unsigned char data[sizeof(T)]; };
    blob b;
    b = *reinterpret_cast<blob*>(&x);
    *reinterpret_cast<blob*>(&x) = *reinterpret_cast<blob*>(&y);
    *reinterpret_cast<blob*>(&y) = b;
}

} // ecell4

#endif /* __ECELL4_UTILS_SWAP_HPP */
