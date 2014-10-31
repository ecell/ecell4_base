#ifndef REACTION_RECORDER_HPP
#define REACTION_RECORDER_HPP

template<typename Trr_>
class ReactionRecorder
{
public:
    typedef Trr_ reaction_record_type;

public:
    virtual ~ReactionRecorder() {}

    virtual void operator()(reaction_record_type const& rec) = 0;
};

#endif /* REACTION_RECORDER_HPP */
