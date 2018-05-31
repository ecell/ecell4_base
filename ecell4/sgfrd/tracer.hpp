#ifndef ECELL4_SGFRD_TRACER
#define ECELL4_SGFRD_TRACER
#include <boost/format.hpp>
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
        std::ofstream f0(fnames[0].c_str()); f0.close();
        std::ofstream f1(fnames[1].c_str()); f1.close();
    }
    basic_tracer(const std::string& fname, const std::size_t indent_size)
        : indent_(0), indent_size_(indent_size), current(0)
    {
        fnames[0] = fname + std::string("01.log");
        fnames[1] = fname + std::string("02.log");
        std::ofstream f0(fnames[0].c_str()); f0.close();
        std::ofstream f1(fnames[1].c_str()); f1.close();
    }
    ~basic_tracer(){}

    void indent()   throw() {indent_ += 1; return;}
    void unindent() throw() {indent_ -= 1; return;}

    void write(const std::string& tr)
    {
        {
            std::ofstream ofs(fnames.at(current).c_str(),
                std::ios_base::in | std::ios_base::out | std::ios_base::ate);
            const std::string idt(indent_size_ * indent_, ' ');
            ofs << idt << tr << std::endl;
            ofs.close();
        }

        std::ifstream ifs(fnames.at(current).c_str());
        ifs.seekg(0, std::ios::beg);
        std::ifstream::streampos init = ifs.tellg();
        ifs.seekg(0, std::ios::end);
        std::ifstream::streampos last = ifs.tellg();
        ifs.close();

        std::ifstream::streampos sz = last - init;
        if(sz > 50000000)
        {
            current = (current == 0) ? 1 : 0;
            // clear the next file
            std::ofstream nxt(fnames.at(current).c_str(), std::ios_base::trunc);
            nxt.close();
        }
        return;
    }

    template<typename formT, typename T1>
    void write(const formT& tr, const T1& a1)
    {
        this->write((boost::format(tr) % a1).str());
    }
    template<typename formT, typename T1, typename T2>
    void write(const formT& tr, const T1& a1, const T2& a2)
    {
        this->write((boost::format(tr) % a1 % a2).str());
    }
    template<typename formT, typename T1, typename T2, typename T3>
    void write(const formT& tr, const T1& a1, const T2& a2, const T3& a3)
    {
        this->write((boost::format(tr) % a1 % a2 % a3).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4>
    void write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4)
    {
        this->write((boost::format(tr) % a1 % a2 % a3 % a4).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5>
    void write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5)
    {
        this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
    void write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6)
    {
        this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5 % a6).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
    void write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6, const T7& a7)
    {
        this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5 % a6 % a7).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
    void write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6, const T7& a7, const T8& a8)
    {
        this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5 % a6 % a7 % a8).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
    void write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6, const T7& a7, const T8& a8, const T9& a9)
    {
        this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5 % a6 % a7 % a8 % a9).str());
    }
    template<typename formT, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
    void write(const formT& tr, const T1& a1, const T2& a2, const T3& a3, const T4& a4, const T5& a5, const T6& a6, const T7& a7, const T8& a8, const T9& a9, const T10& a10)
    {
        this->write((boost::format(tr) % a1 % a2 % a3 % a4 % a5 % a6 % a7 % a8 % a9 % a10).str());
    }

  private:
    const std::size_t indent_size_;
    std::size_t       indent_;
    std::size_t       current;
    boost::array<std::string, 2> fnames;
};

typedef basic_tracer<char, std::char_traits<char> > tracer;

template<typename charT = char, typename traitsT = std::char_traits<charT> >
class scope
{
  public:
    typedef charT   char_type;
    typedef traitsT traits_type;
    typedef basic_tracer<char_type, traits_type> tracer_type;

  public:
    scope(tracer_type& trc)
      : tracer_(trc), name_("")
    {
        tracer_.write("{");
        tracer_.indent();
    }
    scope(tracer_type& trc, const std::string& name)
      : tracer_(trc), name_(name)
    {
        tracer_.write("%s {", name_);
        tracer_.indent();
    }
    ~scope()
    {
        tracer_.unindent();
        tracer_.write("}");
    }

    std::string const& name() const throw() {return name_;}

  private:
    tracer_type& tracer_;
    const std::string name_;
};

typedef scope</*boost::chrono::nanoseconds, */ char, std::char_traits<char> > scope_ns;
typedef scope</*boost::chrono::microseconds,*/ char, std::char_traits<char> > scope_us;
typedef scope</*boost::chrono::milliseconds,*/ char, std::char_traits<char> > scope_ms;
typedef scope</*boost::chrono::seconds,     */ char, std::char_traits<char> > scope_s;
typedef scope</*boost::chrono::minutes,     */ char, std::char_traits<char> > scope_m;
typedef scope</*boost::chrono::hours,       */ char, std::char_traits<char> > scope_h;

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
