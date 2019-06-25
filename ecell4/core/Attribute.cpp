#include "Attribute.hpp"

#include <algorithm>


namespace ecell4
{

Attribute::Attribute()
    : attributes_()
{
    ; // do nothing
}

Attribute::Attribute(const Attribute& another)
    : attributes_(another.attributes_)
{
    ;
}

Attribute::Attribute(const Attribute::container_type& attr)
    : attributes_(attr)
{
    ;
}

Attribute& Attribute::operator=(const Attribute& another)
{
    attributes_ = another.attributes_;
    return *this;
}

std::vector<Attribute::value_type> Attribute::list_attributes() const
{
    std::vector<value_type> retval;
    for (container_type::const_iterator
        i(attributes_.begin()); i != attributes_.end(); ++i)
    {
        retval.push_back(*i);
    }
    return retval;
}

Attribute::mapped_type Attribute::get_attribute(const key_type& key) const
{
    container_type::const_iterator
        i(attributes_.find(key));
    if (i == attributes_.end())
    {
        std::ostringstream message;
        message << "attribute [" << key << "] not found";
        throw NotFound(message.str()); // use boost::format if it's allowed
    }

    return (*i).second;
}

void Attribute::set_attributes(const Attribute& attr)
{
    //XXX: Deprecate me!!!
    attributes_ = attr.attributes_;
}

void Attribute::overwrite_attributes(const Attribute& attr)
{
    const container_type& attrs(attr.attributes_);
    for (container_type::const_iterator i(attrs.begin());
        i != attrs.end(); ++i)
    {
        this->set_attribute((*i).first, (*i).second);
    }
}

void Attribute::remove_attribute(const key_type& key)
{
    container_type::iterator
        i(attributes_.find(key));
    if (i == attributes_.end())
    {
        std::ostringstream message;
        message << "attribute [" << key << "] not found";
        throw NotFound(message.str()); // use boost::format if it's allowed
    }

    attributes_.erase(i);
}

bool Attribute::has_key(const key_type& key) const
{
    return (attributes_.find(key) != attributes_.end());
}

template <>
Real Attribute::get_attribute_as<Real>(const key_type& key) const
{
    mapped_type val = get_attribute(key);
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
Integer Attribute::get_attribute_as<Integer>(const key_type& key) const
{
    mapped_type val = get_attribute(key);
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
void Attribute::set_attribute<const char*>(const key_type& key, const char* value)
{
    set_attribute(key, std::string(value));
}

template <>
void Attribute::set_attribute<Real>(const key_type& key, const Real value)
{
    set_attribute(key, Quantity<Real>(value));
}

template <>
void Attribute::set_attribute<Integer>(const key_type& key, const Integer value)
{
    set_attribute(key, Quantity<Integer>(value));
}

} // ecell4
