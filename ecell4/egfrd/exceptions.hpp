#ifndef ECELL4_EGFRD_EXCEPTIONS_HPP
#define ECELL4_EGFRD_EXCEPTIONS_HPP

#include <ecell4/core/exceptions.hpp>

namespace ecell4
{
namespace egfrd
{

class PropagationError
    : public ecell4::Exception
{
public:

    PropagationError(const std::string& str)
        : str_(str)
    {
        ;
    }

    virtual ~PropagationError() throw()
    {
        ;
    }

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:

    std::string str_;
};

class NoSpace
    : public ecell4::Exception
{
public:

    NoSpace()
        : str_()
    {
        ;
    }

    NoSpace(const std::string& str)
        : str_(str)
    {
        ;
    }

    virtual ~NoSpace() throw()
    {
        ;
    }

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:

    std::string str_;
};

} // egfrd
} // ecell4
#endif /* EXCEPTIONS_HPP */
