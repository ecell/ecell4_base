#ifndef ECELL4_SGFRD_EXPECTED
#define ECELL4_SGFRD_EXPECTED
#include <type_traits>
#include <utility>
#include <stdexcept>
#include <boost/variant.hpp>
#include <boost/optional.hpp>

namespace ecell4
{

template<typename T>
struct success
{
    typedef T value_type;

    success() = default;
    ~success() = default;
    success(const success&) = default;
    success(success&&)      = default;
    success& operator=(const success&) = default;
    success& operator=(success&&)      = default;

    explicit success(const value_type& v)
        noexcept(noexcept(std::is_nothrow_copy_constructible<value_type>::value))
        : value(v)
    {}
    explicit success(value_type&& v)
        noexcept(noexcept(std::is_nothrow_move_constructible<value_type>::value))
        : value(std::move(v))
    {}

    value_type value;
};

template<>
struct success<void>
{
    typedef void value_type;

    success() = default;
    ~success() = default;
    success(const success&) = default;
    success(success&&)      = default;
    success& operator=(const success&) = default;
    success& operator=(success&&)      = default;
};

template<typename T>
struct failure
{
    typedef T value_type;

    failure() = default;
    ~failure() = default;
    failure(const failure&) = default;
    failure(failure&&)      = default;
    failure& operator=(const failure&) = default;
    failure& operator=(failure&&)      = default;

    explicit failure(const value_type& v)
        noexcept(noexcept(std::is_nothrow_copy_constructible<value_type>::value))
        : value(v)
    {}
    explicit failure(value_type&& v)
        noexcept(noexcept(std::is_nothrow_move_constructible<value_type>::value))
        : value(std::move(v))
    {}

    value_type value;
};

template<>
struct failure<void>
{
    typedef void value_type;

    failure() = default;
    ~failure() = default;
    failure(const failure&) = default;
    failure(failure&&)      = default;
    failure& operator=(const failure&) = default;
    failure& operator=(failure&&)      = default;
};

template<typename T>
inline
success<typename std::remove_cv<typename std::remove_reference<T>::type>::type>
ok(T&& v)
{
    return success<
        typename std::remove_cv<typename std::remove_reference<T>::type>::type
        >(std::forward<T>(v));
}
inline success<void> ok()
{
    return success<void>();
}

template<typename T>
inline
failure<typename std::remove_cv<typename std::remove_reference<T>::type>::type>
err(T&& v)
{
    return failure<
        typename std::remove_cv<typename std::remove_reference<T>::type>::type
        >(std::forward<T>(v));
}
inline failure<void> err()
{
    return failure<void>();
}

template<typename T, typename E>
class expected
{
  public:
    typedef T value_type;
    typedef E error_type;
    typedef success<value_type> success_type;
    typedef failure<error_type> failure_type;

  public:

    expected()  = delete; // expected should be initialized
    ~expected() = default;
    expected(const expected& v) = default;
    expected(expected&& v)      = default;
    expected& operator=(const expected& v) = default;
    expected& operator=(expected&& v)      = default;

    expected(const success_type& s) : result_(s) {}
    expected(const failure_type& f) : result_(f) {}
    expected(success_type&& s) : result_(std::move(s)) {}
    expected(failure_type&& f) : result_(std::move(f)) {}
    expected& operator=(const success_type& v){result_ = v; return *this;}
    expected& operator=(const failure_type& e){result_ = e; return *this;}
    expected& operator=(success_type&& v){result_ = std::move(v); return *this;}
    expected& operator=(failure_type&& e){result_ = std::move(e); return *this;}

    const value_type& unwrap() const
    {
        if(!this->is_ok())
        {
            throw std::runtime_error("ecell4::expected::unwrap: not ok");
        }
        return boost::get<success_type>(result_).value;
    }

    const error_type& unwrap_error() const
    {
        if(!this->is_err())
        {
            throw std::runtime_error(
                    "ecell4::expected::unwrap_error: not an error");
        }
        return boost::get<failure_type>(this->result_).value;
    }

    bool is_ok()    const noexcept {return this->result_.which() == 0;}
    bool is_err()   const noexcept {return this->result_.which() == 1;}
    operator bool() const noexcept {return this->is_ok();}

    boost::optional<const value_type&> ok()  const noexcept
    {
        if(this->is_ok())
        {
            return boost::get<success_type>(this->result_).value;
        }
        return boost::none;
    }
    boost::optional<const error_type&> err() const noexcept
    {
        if(this->is_err())
        {
            return boost::get<failure_type>(this->result_).value;
        }
        return boost::none;
    }

  private:

    boost::variant<success_type, failure_type> result_;
};

} // ecell4
#endif// ECELL4_SGFRD_EXPECTED
