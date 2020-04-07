#ifndef ECELL4_EXCEPTIONS_HPP
#define ECELL4_EXCEPTIONS_HPP

#include <exception>
#include <stdexcept>
#include <sstream>
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

namespace detail
{
inline std::string concat_arguments_to_string_impl(std::ostringstream& oss)
{
    return oss.str();
}
template<typename T, typename ... Ts>
std::string concat_arguments_to_string_impl(
        std::ostringstream& oss, T&& head, Ts&& ... tail)
{
    oss << std::forward<T>(head);
    return concat_arguments_to_string_impl(oss, std::forward<Ts>(tail)...);
}
template<typename ... Ts>
std::string concat_arguments_to_string(Ts&& ... args)
{
    std::ostringstream oss;
    return concat_arguments_to_string_impl(oss, std::forward<Ts>(args)...);
}
} // detail

template<class Exception, typename ... Ts>
[[noreturn]] void throw_exception(Ts&& ... args)
{
    throw Exception(detail::concat_arguments_to_string(std::forward<Ts>(args)...));
}

} // ecell4
#endif /* ECELL4_EXCEPTIONS_HPP */
