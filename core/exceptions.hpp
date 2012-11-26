#ifndef __EXCEPTIONS_HPP
#define __EXCEPTIONS_HPP

#include <exception>
#include <stdexcept>


namespace ecell4
{

class Exception
    : public std::exception
{
public:

    Exception()
    {
        ;
    }

    virtual ~Exception() throw()
    {
        ;
    }

    virtual const char* what() const throw()
    {
        return "";
    }
};

class NotFound
    : public Exception
{
public:

    NotFound(std::string const& str)
        : str_(str)
    {
        ;
    }

    virtual ~NotFound() throw()
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

class AlreadyExists
    : public Exception
{
public:

    AlreadyExists(std::string const& str)
        : str_(str)
    {
        ;
    }

    virtual ~AlreadyExists() throw()
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

class NotImplemented
    : public Exception
{
public:

    NotImplemented(std::string const& str)
        : str_(str)
    {
        ;
    }

    virtual ~NotImplemented() throw()
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

class NotSupported
    : public Exception
{
public:

    NotSupported(std::string const& str)
        : str_(str)
    {
        ;
    }

    virtual ~NotSupported() throw()
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

}

#endif /* __EXCEPTIONS_HPP */
