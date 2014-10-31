#ifndef REFERENCE_OR_INSTANCE_HPP
#define REFERENCE_OR_INSTANCE_HPP

template<typename T_>
class reference_or_instance
{
public:
    reference_or_instance(T_& r): type_(REF)
    {
        ref_ = &r;
    }

    reference_or_instance(T_ const& i, int)
        : type_(INSTANCE), ref_(new (instance_) T_(i)) {}

    reference_or_instance(): type_(INSTANCE), ref_(new (instance_) T_()) {}

    reference_or_instance(reference_or_instance const& that): type_(that.type_)
    {
        switch (type_)
        {
        case REF:
            ref_ = that.ref_;
            break;
        case INSTANCE:
            ref_ = new (instance_) T_(static_cast<T_ const&>(that));
            break;
        }
    }

    ~reference_or_instance()
    {
        if (type_ == INSTANCE)
        {
            ref_->~T_();
        }
    }

    operator T_&()
    {
        return *ref_;
    }

    operator T_ const&() const
    {
        return *ref_;
    }

private:
    enum {
        REF,
        INSTANCE
    } type_;

    char instance_[sizeof(T_)];
    T_* ref_;
};

#endif /* REFERENCE_OR_INSTANCE_HPP */
