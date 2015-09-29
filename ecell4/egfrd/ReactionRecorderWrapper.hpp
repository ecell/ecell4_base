#ifndef REACTION_RECORDER_WRAPPER_HPP
#define REACTION_RECORDER_WRAPPER_HPP

#include <boost/shared_ptr.hpp>
#include <ecell4/core/ReactionRule.hpp>
#include <ecell4/core/Identifier.hpp>
#include "ReactionRecorder.hpp"
#include "ReactionRecord.hpp"


template<typename Trr_>
class ReactionRecorderWrapper
    : public ReactionRecorder<Trr_>
{
public:

    typedef ReactionRecorder<Trr_> base_type;

    typedef typename base_type::reaction_record_type reaction_record_type;
    typedef typename reaction_record_type::particle_id_type particle_id_type;
    typedef typename reaction_record_type::reaction_rule_id_type reaction_rule_id_type;
    typedef typename reaction_record_type::reactants_type reactants_type;
    typedef typename reaction_record_type::products_type products_type;

    typedef reaction_record_type reaction_info_type;

public:

    ReactionRecorderWrapper()
        : backend_()
    {
        ;
    }

    virtual ~ReactionRecorderWrapper()
    {
        ;
    }

    virtual void operator()(reaction_record_type const& rec)
    {
        if (backend_)
        {
            (*backend_)(rec);
        }

        // last_reactions_.push_back(rec.reaction_rule_id());
        last_reactions_.push_back(std::make_pair(rec.reaction_rule_id(), rec));
    }

    const std::vector<std::pair<ecell4::ReactionRule, reaction_info_type> >& last_reactions() const
    {
        return last_reactions_;
    }

    void clear()
    {
        last_reactions_.clear();
    }

    boost::shared_ptr<base_type> const& backend() const
    {
        return backend_;
    }

    boost::shared_ptr<base_type>& backend()
    {
        return backend_;
    }

protected:

    std::vector<std::pair<ecell4::ReactionRule, reaction_info_type> > last_reactions_;
    boost::shared_ptr<base_type> backend_;
};


#endif /* REACTION_RECORDER_WRAPPER_HPP */
