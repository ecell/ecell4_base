#ifndef GFRD_POLYGON_CIRCULAR_ITERATION
#define GFRD_POLYGON_CIRCULAR_ITERATION
#include <iterator>
#include <cstddef>
#include <cassert>

namespace detail
{

template<std::size_t I_size, std::size_t I_times>
struct increment_repeat_impl
{
    static std::size_t apply(const std::size_t i); 
};

template<std::size_t I_size, std::size_t I_times>
struct decrement_repeat_impl
{
    static std::size_t apply(const std::size_t i); 
};

}//detail

template<std::size_t I_size>
struct circular_iteration
{
    static std::size_t increment(const std::size_t i)
    {
        return (i < I_size - 1) ? (i + 1) : (i + 1 - I_size);
    }
    static std::size_t decrement(const std::size_t i)
    {
        return (0 < i) ? (i - 1) : (i - 1 + I_size);
    }

    template<std::size_t I_times>
    static std::size_t increment_repeat(const std::size_t i)
    {
        return detail::increment_repeat_impl<I_size, I_times>::apply(i);
    }

    template<std::size_t I_times>
    static std::size_t decrement_repeat(const std::size_t i)
    {
        return detail::decrement_repeat_impl<I_size, I_times>::apply(i);
    }
};

namespace detail
{

template<std::size_t I_size, std::size_t I_times>
inline std::size_t increment_repeat_impl<I_size, I_times>::apply(const std::size_t i)
{
    return increment_repeat_impl<I_size, I_times - 1>::apply(
               circular_iteration<I_size>::increment(i));
}

template<std::size_t I_size, std::size_t I_times>
inline std::size_t decrement_repeat_impl<I_size, I_times>::apply(const std::size_t i)
{
    return decrement_repeat_impl<I_size, I_times - 1>::apply(
               circular_iteration<I_size>::decrement(i));
}

template<std::size_t I_size>
struct increment_repeat_impl<I_size, 0>
{
    static std::size_t apply(const std::size_t i){return i;}
};

template<std::size_t I_size>
struct decrement_repeat_impl<I_size, 0>
{
    static std::size_t apply(const std::size_t i){return i;}
};

}//detail

template<typename T_iter>
class circular_iterator
{
  public:
    typedef T_iter                                      iterator_type;
    typedef std::iterator_traits<iterator_type>         iterator_traits;
    typedef typename iterator_traits::difference_type   difference_type;
    typedef typename iterator_traits::value_type        value_type;
    typedef typename iterator_traits::pointer           pointer;
    typedef typename iterator_traits::reference         reference;
    typedef typename iterator_traits::iterator_category iterator_category;

  public:
    circular_iterator(T_iter begin, T_iter end, T_iter start)
        : round_times(0), begin_(begin), end_(end), start_(start), iter_(start)
    {
        assert(std::distance(begin, start) >= 0);
        assert(std::distance(start, end) >= 0);
    }
    ~circular_iterator(){}

    reference operator*()  const throw() {return *iter_;}
    pointer   operator->() const throw() {return &(*iter_);}

    iterator_type& operator++()    throw();
    iterator_type& operator++(int) throw();
    iterator_type& operator--()    throw();
    iterator_type& operator--(int) throw();

    bool operator==(const iterator_type& i) const throw() {return (iter_ == i);}
    bool operator!=(const iterator_type& i) const throw() {return (iter_ != i);}

    int  round() const throw() {return round_times;}
    bool end()   const throw() {return (round_times == 1) && (iter_ == start_);}

  private:
    int round_times;
    const iterator_type begin_;
    const iterator_type start_;
    const iterator_type end_;
          iterator_type iter_;
};

template<typename T_iter>
typename circular_iterator<T_iter>::iterator_type&
circular_iterator<T_iter>::operator++() throw()
{
    (this->iter_)++;
    if(iter_ == this->end_){iter_ = this->begin_; round_times++;}
    return this->iter_;
}

template<typename T_iter>
typename circular_iterator<T_iter>::iterator_type&
circular_iterator<T_iter>::operator++(int) throw()
{
    ++(this->iter_);
    if(iter_ == this->end_){iter_ = this->begin_; round_times++;}
    return this->iter_;
}

template<typename T_iter>
typename circular_iterator<T_iter>::iterator_type&
circular_iterator<T_iter>::operator--() throw()
{
    if(iter_ == this->begin_){iter_ = this->end_; round_times--;}
    (this->iter_)--;
    return this->iter_;
}

template<typename T_iter>
typename circular_iterator<T_iter>::iterator_type&
circular_iterator<T_iter>::operator--(int) throw()
{
    if(iter_ == this->begin_){iter_ = this->end_; round_times--;}
    --(this->iter_);
    return this->iter_;
}

#endif /* lGFRD_POLYGON_CIRCULAR_ITERATION */
