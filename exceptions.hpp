#ifndef EXCEPTIONS_HPP
#define EXCEPTIONS_HPP

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

#endif /* EXCEPTIONS_HPP */
