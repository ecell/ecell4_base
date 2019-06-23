#ifndef ECELL4_SGFRD_STATISTICS
#define ECELL4_SGFRD_STATISTICS
#include <map>

namespace ecell4
{
namespace sgfrd
{

template<typename Kind>
class statistics
{
  public:
    typedef std::size_t count_type;
    typedef std::map<Kind, count_type> container_type;

  public:

    void add_count(const Kind& k)
    {
        if(this->values_.count(k) == 0)
        {
            this->values_[k] = 0;
        }
        this->values_[k] += 1;
        return;
    }

    std::size_t total() const
    {
        std::size_t tot = 0;
        for(typename container_type::const_iterator
                i(values_.begin()), e(values_.end()); i != e; ++i)
        {
            tot += i->second;
        }
        return tot;
    }

    std::size_t show_count(const Kind& k) const
    {
        if(this->values_.count(k) == 0)
        {
            return 0;
        }
        return values_.find(k)->second;
    }

    double show_percent(const Kind& k) const
    {
        return (100.0 * this->show_count(k)) / this->total();
    }

    std::vector<Kind> list_all_kinds() const
    {
        std::vector<Kind> ks;
        for(typename container_type::const_iterator
                i(values_.begin()), e(values_.end()); i != e; ++i)
        {
            ks.push_back(i->first);
        }
        return ks;
    }

  private:
    container_type values_;
};

enum ReactionKind
{
    SingleFirstOrder,
    SingleFirstOrderFailed,
    PairFirstOrder,
    PairFirstOrderFailed,
    PairSecondOrder,
    PairSecondOrderFailed,
    MultiFirstOrder,
    MultiSecondOrder
};

enum EventFired
{
    FireSingleCircular,
    FireSingleConical,
    FirePair,
    FireMulti,
    FireBirth
};

enum MultiReason
{
    SingleConicalFailed,
    PairFailed
};

#ifndef STATISTICS_OFF
#define STAT(x) x
#else
#define STAT(x) /**/
#endif

} // sgfrd
} // ecell4
#endif // ECELL4_SGFRD_STATISTICS
