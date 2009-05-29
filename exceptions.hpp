#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

#include <string>
#include <exception>

class not_found: public std::exception
{
public:
    not_found(std::string const& str): str_(str) {}

    virtual ~not_found() throw() {}

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:
    std::string str_;
};

class already_exists: public std::exception
{
public:
    already_exists(std::string const& str): str_(str) {}

    virtual ~already_exists() throw() {}

    virtual const char* what() const throw()
    {
        return str_.c_str();
    }

private:
    std::string str_;
};

#endif /* EXCEPTIONS_HPP */
