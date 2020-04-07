#include "Attribute.hpp"

#include <algorithm>


namespace ecell4
{

Attribute::Attribute()
    : data_()
{
    ; // do nothing
}

Attribute::Attribute(const Attribute& another)
    : data_(another.data_)
{
    ;
}

Attribute::Attribute(const Attribute::container_type& attr)
    : data_(attr)
{
    ;
}

Attribute& Attribute::operator=(const Attribute& another)
{
    data_ = another.data_;
    return *this;
}

std::vector<Attribute::value_type> Attribute::values() const
{
    std::vector<value_type> retval;
    for (container_type::const_iterator
        i(data_.begin()); i != data_.end(); ++i)
    {
        retval.push_back(*i);
    }
    return retval;
}

Attribute::mapped_type Attribute::get(const key_type& key) const
{
    container_type::const_iterator
        i(data_.find(key));
    if (i == data_.end())
    {
        throw_exception<NotFound>("attribute [", key, "] not found");
    }

    return (*i).second;
}

void Attribute::clear()
{
    data_.clear();
}

void Attribute::overwrite(const Attribute& attr)
{
    const container_type& attrs(attr.data_);
    for (container_type::const_iterator i(attrs.begin());
        i != attrs.end(); ++i)
    {
        this->set((*i).first, (*i).second);
    }
}

void Attribute::remove(const key_type& key)
{
    container_type::iterator
        i(data_.find(key));
    if (i == data_.end())
    {
        throw_exception<NotFound>("attribute [", key, "] not found");
    }

    data_.erase(i);
}

bool Attribute::has_key(const key_type& key) const
{
    return (data_.find(key) != data_.end());
}

template <>
Real Attribute::get_as<Real>(const key_type& key) const
{
    mapped_type val = get(key);
    if (Quantity<Real>* x = boost::get<Quantity<Real> >(&val))
    {
        return (*x).magnitude;
    }
    else if (Quantity<Integer>* x = boost::get<Quantity<Integer> >(&val))
    {
        return static_cast<Real>((*x).magnitude);
    }
    else if (std::string* x = boost::get<std::string>(&val))
    {
        return std::atof((*x).c_str());
    }
    throw NotSupported("An attribute has incorrect type. Real is expected");
}

template <>
Integer Attribute::get_as<Integer>(const key_type& key) const
{
    mapped_type val = get(key);
    if (Quantity<Integer>* x = boost::get<Quantity<Integer> >(&val))
    {
        return (*x).magnitude;
    }
    else if (std::string* x = boost::get<std::string>(&val))
    {
        return std::atoi((*x).c_str());
    }
    throw NotSupported("An attribute has incorrect type. Integer is expected");
}

template <>
void Attribute::set<const char*>(const key_type& key, const char* value)
{
    set(key, std::string(value));
}

template <>
void Attribute::set<Real>(const key_type& key, const Real value)
{
    set(key, Quantity<Real>(value));
}

template <>
void Attribute::set<Integer>(const key_type& key, const Integer value)
{
    set(key, Quantity<Integer>(value));
}

} // ecell4
