#ifndef ECELL4_EXCEPTIONS_HPP
#define ECELL4_EXCEPTIONS_HPP

#include <exception>
#include <stdexcept>

#include <string>

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

    NotFound(const std::string& str)
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

    AlreadyExists(const std::string& str)
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

    NotImplemented(const std::string& str)
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

    NotSupported(const std::string& str)
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

class IllegalState
    : public Exception
{
public:

    IllegalState(const std::string& str)
        : str_(str)
    {
        ;
    }

    virtual ~IllegalState() throw()
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

class IllegalArgument
    : public Exception
{
public:

    IllegalArgument(const std::string& str)
        : str_(str)
    {
        ;
    }


    virtual ~IllegalArgument() throw()
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

#endif /* ECELL4_EXCEPTIONS_HPP */
