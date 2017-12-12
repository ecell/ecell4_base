#ifndef ECELL4_SGFRD_EXPECTED
#define ECELL4_SGFRD_EXPECTED
#include <boost/variant.hpp>

namespace ecell4
{

template<typename T, typename E>
class expected
{
  public:
    typedef T value_type;
    typedef E error_type;

  public:

    explicit expected(const value_type& v) : result_(v) {}
    explicit expected(const error_type& e) : result_(e) {}
    expected& operator=(const value_type& v){this->result_ = v; return *this;}
    expected& operator=(const error_type& e){this->result_ = e; return *this;}

    expected(const expected& v) : result_(v.result_) {}
    expected& operator=(const expected& v)
    {
        this->result_ = v.result_;
        return *this;
    }

    const value_type& unwrap() const throw()
    {
        assert(this->is_ok());
        return boost::get<value_type>(result_);
    }
    const value_type& unwrap_or(const value_type& opt) const throw()
    {
        if(this->is_ok()) {return boost::get<value_type>(this->result_);}
        else {return opt;}
    }

    const error_type& unwrap_error() const throw()
    {
        assert(this->is_err());
        return boost::get<error_type>(this->result_);
    }
    const error_type& unwrap_error_or(const error_type& opt) const throw()
    {
        if(this->is_err()) {return boost::get<error_type>(this->result_);}
        else {return opt;}
    }

    bool is_ok()    const throw() {return this->result_.which() == 0;}
    bool is_err()   const throw() {return this->result_.which() == 1;}
    operator bool() const throw() {return this->is_ok();}

    boost::optional<const value_type&> ok() const throw()
    {
        if(this->is_ok()) {return boost::get<value_type>(this->result_);}
        else {return boost::none;}
    }
    boost::optional<const error_type&> err() const throw()
    {
        if(this->is_err()) {return boost::get<error_type>(this->result_);}
        else {return boost::none;}
    }

  private:

    boost::variant<T, E> result_;
};

} // ecell4
#endif// ECELL4_SGFRD_EXPECTED
