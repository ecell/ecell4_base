#ifndef REACTION_RECORDER_HPP
#define REACTION_RECORDER_HPP

namespace ecell4
{
namespace egfrd
{
template<typename Trr_>
class ReactionRecorder
{
public:
    typedef Trr_ reaction_record_type;

public:
    virtual ~ReactionRecorder() {}

    virtual void operator()(reaction_record_type const& rec) = 0;
};

} // egfrd
} // ecell4
#endif /* REACTION_RECORDER_HPP */
