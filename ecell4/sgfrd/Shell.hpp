#ifndef ECELL4_SGFRD_SHELL
#define ECELL4_SGFRD_SHELL
#include <ecell4/core/Real3.hpp>
#include <boost/shared_ptr.hpp>

namespace ecell4
{
namespace sgfrd
{

template<typename T_shape>
class Shell
{
  public:
    typedef T_shape shape_type;

  public:
    Shell(){}
    ~Shell(){}

    Shell(const domain_id_type& domain_id, const shape_type& shape)
        : domain_id_(domain_id), shape_(shape)
    {}

    Real3 const& position() const {return shape_.position();}
    Real3      & position()       {return shape_.position();}

    Real  size() const {return shape_.size();}
    Real& size()       {return shape_.size();}

    DomainID      & domain_id()       {return domain_id_;}
    DomainID const& domain_id() const {return domain_id_;}

    shape_type      & shape()       {return shape_;}
    shape_type const& shape() const {return shape_;}

  private:

    domain_id_type domain_id_;
    shape_type     shape_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SHELL */
