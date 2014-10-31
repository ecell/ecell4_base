#ifndef OBJECTMATRIX_PEER_RANGE_CONVERTERS_HPP
#define OBJECTMATRIX_PEER_RANGE_CONVERTERS_HPP

#include <boost/python.hpp>
#include "generator.hpp"

namespace peer { namespace wrappers {

template<typename Tvalue_>
class pyiterator_generator: public abstract_limited_generator<Tvalue_>
{
public:
    typedef Tvalue_ value_type;

public:
    pyiterator_generator(boost::python::handle<> iter)
        : iter_(iter), advanced_(false) {}

    virtual ~pyiterator_generator() {}

    virtual bool valid() const
    {
        const_cast<pyiterator_generator*>(this)->fetch();
        return last_;
    }

    virtual value_type operator()()
    {
        fetch();
        if (!last_)
        {
            return value_type();
        }
        advanced_ = false;
        return boost::python::extract<value_type>(last_.get())();
    }

    bool operator==(pyiterator_generator const& rhs) const
    {
        return (!last_ && !rhs.last_) || iter_ == rhs.iter_;
    }

private:
    void fetch()
    {
        if (iter_ && !advanced_)
        {
            last_ = boost::python::handle<>(
                    boost::python::allow_null(
                        PyIter_Next(iter_.get())));
            if (!last_)
            {
                iter_.reset();
            }
            advanced_ = true;
        }
    }

protected:
    boost::python::handle<> iter_; 
    bool advanced_;
    boost::python::handle<> last_;
};

} } // namespace peer::wrappers

template<typename Tvalue>
inline bool valid(peer::wrappers::pyiterator_generator<Tvalue> const& gen)
{
    return gen.valid();
}

#endif /* OBJECTMATRIX_PEER_RANGE_CONVERTERS_HPP */
