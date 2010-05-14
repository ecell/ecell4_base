#ifndef PEER_WRAPPERS_ITERATOR_PYSEQ_ITERATOR_HPP
#define PEER_WRAPPERS_ITERATOR_PYSEQ_ITERATOR_HPP

#include <boost/python.hpp>
#include <boost/iterator/iterator_facade.hpp>

template<typename Tvalue_>
class pyseq_iterator
    : public boost::iterator_facade<
        pyseq_iterator<Tvalue_>, Tvalue_,
        boost::random_access_traversal_tag,
        Tvalue_,
        Py_ssize_t>
{
public:
    typedef Py_ssize_t difference_type;
    typedef Tvalue_ reference;

public:
    pyseq_iterator(boost::python::object seq, difference_type idx = 0)
        : seq_(seq), idx_(idx) {}

    reference dereference() const
    {
        return boost::python::extract<Tvalue_>(
            boost::python::handle<>(
                PySequence_GetItem(seq_.ptr(), idx_)).get())();
    }

    bool equal(pyseq_iterator const& rhs) const
    {
        return seq_ == rhs.seq_ && idx_ == rhs.idx_;
    }

    void increment()
    {
        ++idx_;
    }

    void decrement()
    {
        --idx_;
    }

    void advance(difference_type n)
    {
        idx_ += n; 
    }

    difference_type distance_to(pyseq_iterator const& rhs) const
    {
        return rhs.idx_ - idx_;
    }

protected:
    boost::python::object seq_;
    difference_type idx_;
};

#endif /* PEER_WRAPPERS_ITERATOR_PYSEQ_ITERATOR_HPP */
