#ifndef BASIC_REACTION_RECORDER_IMPL_HPP
#define BASIC_REACTION_RECORDER_IMPL_HPP

#include <vector>
#include "ReactionRecorder.hpp"

template<typename Trr_>
class BasicReactionRecorderImpl: public ReactionRecorder<Trr_>
{
public:
    typedef Trr_ reaction_record_type;

public:
    virtual void operator()(reaction_record_type const& rec)
    {
        records_.push_back(rec);
    }


protected:
    std::vector<reaction_record_type> records_;
};

#endif /* BASIC_REACTION_RECORDER_IMPL_HPP */
