#ifndef UTILS_SWAP_HPP
#define UTILS_SWAP_HPP

#include <cstring>

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

#endif /* UTILS_SWAP_HPP */
