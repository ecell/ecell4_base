#ifndef CONSOLE_LOGGER_HPP
#define CONSOLE_LOGGER_HPP

#include <string>

#include "Logger.hpp"

namespace ecell4
{
namespace egfrd
{

class ConsoleAppender: public LogAppender
{
public:
    typedef LogAppender base_type;

public:
    virtual ~ConsoleAppender();

    virtual void flush();

    virtual void operator()(enum Logger::level lv, char const* name, char const** chunks);
};

} // egfrd
} // ecell4
#endif /* CONSOLE_LOGGER_HPP */
