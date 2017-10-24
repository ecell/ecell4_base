#ifndef ECELL4_SGFRD_TRACER
#define ECELL4_SGFRD_TRACER
#include <boost/format.hpp>
#include <boost/chrono.hpp>
#include <boost/array.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <string>
#include <fstream>

namespace ecell4
{
namespace sgfrd
{

template<typename charT, typename traitsT = std::char_traits<charT> >
class basic_tracer;

template<>
class basic_tracer<char, std::char_traits<char> >
{
  public:
    typedef char char_type;
    typedef std::char_traits<char> traits_type;

  public:
    basic_tracer(const std::string& fname)
        : indent_(0), indent_size_(2), current(0)
    {
        fnames[0] = fname + std::string("01.log");
        fnames[1] = fname + std::string("02.log");
        std::ofstream f0(fnames[0]); f0.close();
        std::ofstream f1(fnames[1]); f1.close();
    }
    basic_tracer(const std::string& fname, const std::size_t indent_size)
        : indent_(0), indent_size_(indent_size), current(0)
    {
        fnames[0] = fname + std::string("01.log");
        fnames[1] = fname + std::string("02.log");
        std::ofstream f0(fnames[0]); f0.close();
        std::ofstream f1(fnames[1]); f1.close();
    }
    ~basic_tracer(){}

    basic_tracer& indent()   throw() {indent_ += 1;}
    basic_tracer& unindent() throw() {indent_ -= 1;}

    basic_tracer& write(const std::string& tr)
    {
        std::ofstream ofs(fnames[current].c_str(),
                std::ios_base::in | std::ios_base::out | std::ios_base::ate);
        const std::string idt(indent_size_ * indent_, ' ');
        ofs << idt << tr << std::endl;

        ofs.seekp(0, std::ios::beg);
        std::ofstream::streampos init = ofs.tellp();
        ofs.seekp(0, std::ios::end);
        std::ofstream::streampos last = ofs.tellp();
        ofs.close();

        std::ofstream::streampos sz = last - init;
        if(sz > 50000000)
        {
            current = (current == 0) ? 1 : 0;
            // clear the next file
            std::ofstream nxt(fnames[current].c_str(), std::ios_base::trunc);
            nxt.close();
        }
        return *this;
    }

    template<typename formT, typename T1>
    basic_tracer& write(const formT& tr, const T1& a1)
    {
        return this->write((boost::format(tr) % a1).str());
    }
    template<typename formT, typename T1, typename T2>
    basic_tracer& write(const formT& tr, const T1& a1, const T2& a2)
    {
        return this->write((boost::format(tr) % a1 % a2).str());
    }
    template<typename formT, typename T1, typename T2, typename T3>
    basic_tracer& write(const formT& tr, const T1& a1, const T2& a2, const T3& a3)
    {
        return this->write((boost::format(tr) % a1 % a2 % a3).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4>
    basic_tracer& write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4)
    {
        return this->write((boost::format(tr) % a1 % a2 % a3 % a4).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5>
    basic_tracer& write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5)
    {
        return this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
    basic_tracer& write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6)
    {
        return this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5 % a6).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
    basic_tracer& write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6, const T7& a7)
    {
        return this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5 % a6 % a7).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
    basic_tracer& write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6, const T7& a7, const T8& a8)
    {
        return this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5 % a6 % a7 % a8).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
    basic_tracer& write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6, const T7& a7, const T8& a8, const T9& a9)
    {
        return this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5 % a6 % a7 % a8 % a9).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
    basic_tracer& write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6, const T7& a7, const T8& a8, const T9& a9, const T10& a10)
    {
        return this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5 % a6 % a7 % a8 % a9 % a10).str());
    }

  private:
    const std::size_t indent_size_;
    std::size_t       indent_;
    std::size_t       current;
    boost::array<std::string, 2> fnames;
};

typedef basic_tracer<char, std::char_traits<char> > tracer;

/* ---------------------------- tracer impl end ---------------------------- */


namespace detail
{

template<typename durT>
struct stringize_impl;

template<>
struct stringize_impl<boost::chrono::nanoseconds>
{
    static const char* invoke() throw() {return "ns";}
};

template<>
struct stringize_impl<boost::chrono::microseconds>
{
    static const char* invoke() throw() {return "us";}
};

template<>
struct stringize_impl<boost::chrono::milliseconds>
{
    static const char* invoke() throw() {return "ms";}
};

template<>
struct stringize_impl<boost::chrono::seconds>
{
    static const char* invoke() throw() {return "s";}
};

template<>
struct stringize_impl<boost::chrono::minutes>
{
    static const char* invoke() throw() {return "m";}
};

template<>
struct stringize_impl<boost::chrono::hours>
{
    static const char* invoke() throw() {return "h";}
};

} // detail

template<typename durT>
inline const char* stringize() throw()
{
    return detail::stringize_impl<durT>::invoke();
}

template<typename durationT, typename charT = char, typename traitsT = std::char_traits<charT> >
class scope
{
  public:
    typedef charT     char_type;
    typedef traitsT   traits_type;
    typedef durationT duration_type;
    typedef basic_tracer<char_type, traits_type> tracer_type;

  public:
    scope(tracer_type& trc)
      : tracer_(trc), start_(boost::chrono::high_resolution_clock::now()),
        name_("")
    {
        tracer_.write("{");
        tracer_.indent();
    }
    scope(tracer_type& trc, const std::string& name)
      : tracer_(trc), start_(boost::chrono::high_resolution_clock::now()),
        name_(name)
    {
        tracer_.write("%s {", name_);
        tracer_.indent();
    }
    ~scope()
    {
        tracer_.unindent();
        tracer_.write("} %1% [%2%]",
            boost::chrono::duration_cast<durationT>(
                boost::chrono::high_resolution_clock::now() - start_).count(),
            stringize<durationT>());
    }

    std::string const& name() const throw() {return name_;}

  private:
    tracer_type& tracer_;
    boost::chrono::steady_clock::time_point start_;
    const std::string name_;
};

typedef scope<boost::chrono::nanoseconds,  char, std::char_traits<char> > scope_ns;
typedef scope<boost::chrono::microseconds, char, std::char_traits<char> > scope_us;
typedef scope<boost::chrono::milliseconds, char, std::char_traits<char> > scope_ms;
typedef scope<boost::chrono::seconds,      char, std::char_traits<char> > scope_s;
typedef scope<boost::chrono::minutes,      char, std::char_traits<char> > scope_m;
typedef scope<boost::chrono::hours,        char, std::char_traits<char> > scope_h;

/* ----------------------------- scope impl end ----------------------------- */

#ifndef ECELL4_SGFRD_NO_TRACE
#define SGFRD_TRACE(x)               x;
#define SGFRD_SCOPE(time, name, trc) BOOST_PP_CAT(scope_, time) BOOST_PP_CAT(scope_, name)(trc, BOOST_PP_STRINGIZE(name));
#else //ECELL4_SGFRD_NO_TRACE
#define SGFRD_TRACE(x)               /**/
#define SGFRD_SCOPE(time, name, trc) /**/
#endif//ECELL4_SGFRD_NO_TRACE

} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_TRACER
