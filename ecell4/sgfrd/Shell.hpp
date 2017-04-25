#ifndef ECELL4_SGFRD_SHELL
#define ECELL4_SGFRD_SHELL
#include <ecell4/core/Real3.hpp>
#include <ecell4/sgfrd/SGFRDEvent.hpp>
#include <boost/shared_ptr.hpp>

namespace ecell4
{
namespace sgfrd
{

template<typename T_shape, typename T_structure_id>
class Shell
{
  public:
    typedef T_shape shape_type;
    typedef T_structure_id structure_id_type;
    typedef DomainID domain_id_type;

  public:
    Shell(){}
    ~Shell(){}

    Shell(const domain_id_type& domain_id, const shape_type& shape,
          const structure_id_type& sid)
        : domain_id_(domain_id), structure_id_(sid), shape_(shape)
    {}

    Real3 const& position() const {return shape_.position();}
    Real3      & position()       {return shape_.position();}

    Real  size() const {return shape_.size();}
    Real& size()       {return shape_.size();}

    domain_id_type      & domain_id()       {return domain_id_;}
    domain_id_type const& domain_id() const {return domain_id_;}

    structure_id_type      & structure_id()       {return structure_id_;}
    structure_id_type const& structure_id() const {return structure_id_;}
    shape_type      & shape()       {return shape_;}
    shape_type const& shape() const {return shape_;}

    std::pair<Real3, structure_id_type> get_surface_position() const
    {
        return std::make_pair(shape_.position(), structure_id_);
    }

  private:

    domain_id_type    domain_id_;
    structure_id_type structure_id_;
    shape_type     shape_;
};

} // sgfrd
} // ecell4
#endif /* ECELL4_SGFRD_SHELL */
