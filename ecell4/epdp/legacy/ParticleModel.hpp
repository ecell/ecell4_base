#ifndef PARTICLE_MODEL_HPP
#define PARTICLE_MODEL_HPP

#include "Defs.hpp"
#include "Model.hpp"
#include "StructureType.hpp"

class ParticleModel: public Model
{
public:
    typedef Model base_type;
    typedef StructureType structure_type_type;
    typedef structure_type_type::identifier_type structure_id_type;

    typedef std::map<structure_id_type, boost::shared_ptr<structure_type_type> > structure_type_map_type;

    typedef select_second<structure_type_map_type::value_type> structure_second_selector_type;

public:
    typedef boost::transform_iterator<structure_second_selector_type,
            structure_type_map_type::const_iterator> structure_type_iterator;
    typedef boost::iterator_range<structure_type_iterator> structure_type_range;

public:
    ParticleModel();

    virtual ~ParticleModel();

    boost::shared_ptr<structure_type_type> get_structure_type_by_id(structure_id_type const& id) const;

    void add_structure_type(boost::shared_ptr<structure_type_type> const& structure);

    structure_type_range get_structure_types() const;

public:
    structure_type_map_type structure_type_map_;
};


#endif /* MODEL_HPP */
