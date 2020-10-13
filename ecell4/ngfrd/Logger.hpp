#ifndef ECELL4_NGFRD_LOGGER_HPP
#define ECELL4_NGFRD_LOGGER_HPP

#include <ecell4/core/type_name_of.hpp>
#include <boost/circular_buffer.hpp>

// to output those
#include <boost/container/small_vector.hpp>

#include <unordered_map>
#include <ostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace ecell4
{
namespace ngfrd
{

namespace logger_detail
{
template<typename T, typename Alloc>
std::ostream& operator<<(std::ostream& os, const std::vector<T, Alloc>& v);
template<typename T, std::size_t N, typename Alloc, typename Options>
std::ostream& operator<<(std::ostream& os, const boost::container::small_vector<T, N, Alloc, Options>& v);
template<typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& p);

template<typename T, typename Alloc>
std::ostream& operator<<(std::ostream& os, const std::vector<T, Alloc>& v)
{
    os << "\"vector<" << utils::type_name_of<T>::value() << ">\":[";
    bool is_front = true;
    for(const auto& elem : v)
    {
        if(!is_front) {os << ", ";} os << v;
        is_front = false;
    }
    os << "]";
    return os;
}
template<typename T, std::size_t N, typename Alloc, typename Options>
std::ostream& operator<<(std::ostream& os,
           const boost::container::small_vector<T, N, Alloc, Options>& v)
{
    os << "\"boost::small_vector<" << utils::type_name_of<T>::value() << ", " << N << ">\":[";
    for(std::size_t i=0; i<v.size(); ++i)
    {
        if(i != 0) {os << ", ";}
        os << v.at(i);
    }
    os << "]";
    return os;
}

template<typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& p)
{
    os << "\"std::pair<" << utils::type_name_of<T>::value() << ", "
       << utils::type_name_of<U>::value() << ">\":{";
    os << p.first << ", " << p.second;
    os << "}";
    return os;
}

} // logger_detail

class Logger
{
public:
    Logger() = default;
    Logger(std::string fname, const std::size_t size)
        : count_(0), indent_(0), buffer_(size), filename_(std::move(fname))
    {}
    Logger(const Logger&)            = delete;
    Logger& operator=(const Logger&) = delete;
    Logger(Logger&&)            = default;
    Logger& operator=(Logger&&) = default;
    ~Logger()
    {
        if(!filename_.empty())
        {
            std::ofstream ofs(this->filename_);
            for(const auto& elem : buffer_)
            {
                ofs << elem << std::endl;
            }
            ofs.close();
        }
    }

    void indent()   {this->indent_ += 2;}
    void unindent() {this->indent_ -= 2;}

    template<typename ... Ts>
    void log(Ts&& ... args)
    {
        this->count_ += 1;

        std::ostringstream oss;
        log_impl(oss, std::string(indent_, ' '), std::forward<Ts>(args)...);
        this->buffer_.push_back(oss.str());

        if(this->count_ == buffer_.capacity())
        {
            this->dump();
            this->count_ = 0;
        }
    }

    void dump()
    {
        std::ofstream ofs(this->filename_);
        for(const auto& elem : buffer_)
        {
            ofs << elem << std::endl;
        }
        ofs.close();
    }

private:

    template<typename Head, typename ... Tail>
    static void log_impl(std::ostream& os, Head&& head, Tail&& ... tail)
    {
        using namespace logger_detail; // to output containers
        os << head;
        return log_impl(os, std::forward<Tail>(tail)...);
    }
    static void log_impl(std::ostream&) {return;}

private:
    std::size_t count_;
    std::size_t indent_;
    boost::circular_buffer<std::string> buffer_;
    std::string filename_;
};

class LoggerManager
{
private:
    static std::unordered_map<std::string, Logger> loggers_;

public:

    static Logger& get_logger()
    {
        if(loggers_.count("default") == 0)
        {
            loggers_.emplace("default", Logger("ngfrd.log", 10000));
        }
        return loggers_.at("default");
    }

    // ------------------------------------------------------------------------

    static void set_logger(const std::string& name, std::string fname, std::size_t N)
    {
        loggers_[name] = Logger(std::move(fname), N);
        return;
    }
    static Logger& get_logger(const std::string& name)
    {
        return loggers_.at(name);
    }
};

class Scope
{
public:
    Scope() = delete;
    Scope(const Scope&) = delete;
    Scope& operator=(const Scope&) = delete;

    Scope(Logger& logger, std::string name, std::string loc)
        : logger_(logger), name_(std::move(name)),
          loc_(loc.substr(loc.rfind('/')+1))
    {
        // remove template typenames (e.g. [with T = int]) from the name
        // because it often become too long to read
        const auto offset = name_.find('[');
        if(offset != std::string::npos)
        {
            name_.erase(name_.begin() + offset, name_.end());
        }
        logger_.log("\"", this->name_, ":", this->loc_, "\":{");
        logger_.indent();
    }
    ~Scope()
    {
        logger_.unindent();
        logger_.log("}");
    }

private:
    std::string name_, loc_;
    Logger& logger_;
};

} // ngfrd
} // ecell4

#if defined(__GNUC__)
#  define ECELL4_NGFRD_LOG_FUNCTION_NAME __PRETTY_FUNCTION__
#elif defined(_MSC_VER)
#  define ECELL4_NGFRD_LOG_FUNCTION_NAME __FUNCSIG__
#else
#  define ECELL4_NGFRD_LOG_FUNCTION_NAME __func__
#endif

#ifndef ECELL4_NGFRD_STRINGIZE
#  define ECELL4_NGFRD_STRINGIZE_AUX(x) #x
#  define ECELL4_NGFRD_STRINGIZE(x)     ECELL4_NGFRD_STRINGIZE_AUX(x)
#endif

#ifndef ECELL4_NGFRD_CONCAT
#define ECELL4_NGFRD_CONCAT_AUX(x, y) x ## y
#define ECELL4_NGFRD_CONCAT(x, y) ECELL4_NGFRD_CONCAT_AUX(x, y)
#endif

#ifdef ECELL4_NGFRD_LOG_DEBUG
#  ifndef ECELL4_NGFRD_LOG_FUNCTION
#    define ECELL4_NGFRD_LOG_FUNCTION()\
         auto& l_o_g_g_e_r_ = ::ecell4::ngfrd::LoggerManager::get_logger();\
         ::ecell4::ngfrd::Scope ECELL4_NGFRD_CONCAT(s_c_o_p_e_, __LINE__) (l_o_g_g_e_r_, ECELL4_NGFRD_LOG_FUNCTION_NAME, __FILE__ ":" ECELL4_NGFRD_STRINGIZE(__LINE__))\
         /**/
#  endif
#  ifndef ECELL4_NGFRD_LOG
#    define ECELL4_NGFRD_LOG(...) l_o_g_g_e_r_.log(__VA_ARGS__)
#  endif
#else // no LOG_DEBUG
#  ifndef ECELL4_NGFRD_LOG_FUNCTION
#    define ECELL4_NGFRD_LOG_FUNCTION() /**/
#  endif
#  ifndef ECELL4_NGFRD_LOG
#    define ECELL4_NGFRD_LOG(...) /**/
#  endif
#endif

#endif
