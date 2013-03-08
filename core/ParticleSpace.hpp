#ifndef __PARTICLE_SPACE_HPP
#define __PARTICLE_SPACE_HPP

#include <cmath>
#include <string.h>
#include <sstream>
#include <cstdio>
#include <cstring>

// #include <gsl/gsl_pow_int.h>

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "functions.hpp"
#include "exceptions.hpp"
#include "Position3.hpp"
#include "Particle.hpp"
#include "Species.hpp"
#include "Space.hpp"

#include <boost/scoped_ptr.hpp>
#include <boost/scoped_array.hpp>

#ifndef H5_NO_NAMESPACE
#ifndef H5_NO_STD
    using std::cout;
    using std::endl;
#endif  // H5_NO_STD
#endif

#include "H5Cpp.h"

#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

const H5std_string DATASET_NAME( "TimePoint" );
const H5std_string MEMBER1( "particle_id" );
const H5std_string MEMBER2( "species_name" );
const H5std_string MEMBER3( "position_x" );
const H5std_string MEMBER4( "position_y" );
const H5std_string MEMBER5( "position_z" );

namespace ecell4
{

Real pow_2(Real const& a);

class ParticleSpace
    : public Space
{
public:

    typedef std::vector<std::pair<ParticleID, Particle> >
    particle_container_type;

    virtual Position3 const& edge_lengths() const
    {
        throw NotImplemented("edge_lengths() not implemented");
    }

    virtual Integer num_particles() const
    {
        throw NotImplemented("num_particles() not implemented");
    }

    virtual Integer num_particles(Species const& species) const
    {
        throw NotImplemented("num_particles() not implemented");
    }

    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles(Species const& species) const
    {
        throw NotImplemented("list_particles() not implemented");
    }

    virtual bool has_particle(ParticleID const& pid) const = 0;
    virtual bool update_particle(ParticleID const& pid, Particle const& p) = 0;
    virtual void remove_particle(ParticleID const& pid) = 0;

    virtual particle_container_type const& particles() const = 0;

    virtual std::pair<ParticleID, Particle>
    get_particle(ParticleID const& pid) const = 0;
    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles() const = 0;

    virtual void save_positions(H5::H5File *file_, Real const& t) = 0;

    Position3 periodic_transpose(
        Position3 const& pos1, Position3 const& pos2) const
    {
        Position3 retval(pos1);
        Position3 const& edges(edge_lengths());
        for (Position3::size_type dim(0); dim < 3; ++dim)
        {
            const Real edge_length(edges[dim]);
            const Real diff(pos2[dim] - pos1[dim]), half(edge_length * 0.5);

            if (diff > half)
            {
                retval[dim] += edge_length;
            }
            else if (diff < -half)
            {
                retval[dim] -= edge_length;
            }
        }
        return retval;
    }

    inline Position3 apply_boundary(Position3 const& pos) const
    {
        return modulo(pos, edge_lengths());
    }

    Real distance_sq(
        Position3 const& pos1, Position3 const& pos2) const
    {
        Real retval(0);
        Position3 const& edges(edge_lengths());
        for (Position3::size_type dim(0); dim < 3; ++dim)
        {
            const Real edge_length(edges[dim]);
            const Real diff(pos2[dim] - pos1[dim]), half(edge_length * 0.5);

            if (diff > half)
            {
                retval += pow_2(diff - edge_length);
            }
            else if (diff < -half)
            {
                retval += pow_2(diff + edge_length);
            }
            else
            {
                retval += pow_2(diff);
            }
        }
        return retval;
    }

    inline Real distance(Position3 const& pos1, Position3 const& pos2) const
    {
        return std::sqrt(distance_sq(pos1, pos2));
    }

    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius) const = 0;
    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius,
        ParticleID const& ignore) const = 0;
    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius,
        ParticleID const& ignore1, ParticleID const& ignore2) const = 0;
};

class ParticleSpaceVectorImpl
    : public ParticleSpace
{
public:

    typedef ParticleSpace::particle_container_type particle_container_type;

    ParticleSpaceVectorImpl(Position3 const& edge_lengths)
    {
        set_edge_lengths(edge_lengths);
    }

    Position3 const& edge_lengths() const
    {
        return edge_lengths_;
    }

    Integer num_particles() const;
    Integer num_particles(Species const& species) const;

    bool has_particle(ParticleID const& pid) const;
    bool update_particle(ParticleID const& pid, Particle const& p);
    void remove_particle(ParticleID const& pid);

    particle_container_type const& particles() const
    {
        return particles_;
    }

    void save_positions(H5::H5File *file_, Real const& t)
    {

    	// Define data structure

    	typedef struct h5_partcle_struct {
    		int h5_particle_id;
    		char h5_species_name[32];
    		double h5_particle_position_x;
    		double h5_particle_position_y;
    		double h5_particle_position_z;
    	} h5_particle_struct;

    	particle_container_type::size_type const np(particles_.size());
        std::cout << np << std::endl;

    	h5_particle_struct h5_p[np];

        // construct data set.
        //boost::scoped_array<h5_particle_struct> particle_id_table(new h5_particle_struct[ particles_.size() ]);

    	for (int i=0; i<np; i++){
    		h5_p[i].h5_particle_id = particles_[i].first;
    		strcpy( h5_p[i].h5_species_name, particles_[i].second.species().name().c_str() );
    		h5_p[i].h5_particle_position_x = particles_[i].second.position()[0];
    		h5_p[i].h5_particle_position_y = particles_[i].second.position()[1];
    		h5_p[i].h5_particle_position_z = particles_[i].second.position()[2];
    	}

//        for (unsigned int i(0); i < particles_.size(); i++)
//        {
//        	particle_id_table[i].h5_particle_id = particles_[i].first;
//        	std::strcpy(particle_id_table[i].h5_species_name, particles_[i].second.species().name().c_str());
//        	particle_id_table[i].h5_particle_position_x = particles_[i].second.position()[0];
//        	particle_id_table[i].h5_particle_position_y = particles_[i].second.position()[1];
//        	particle_id_table[i].h5_particle_position_z = particles_[i].second.position()[2];
//        }


    	//H5::Exception::dontPrint();

    	CompType mtype( sizeof(h5_particle_struct) );
        mtype.insertMember( MEMBER1, HOFFSET(h5_particle_struct, h5_particle_id), PredType::NATIVE_INT);
        mtype.insertMember( MEMBER2, HOFFSET(h5_particle_struct, h5_species_name), StrType(PredType::C_S1, 32) );
        mtype.insertMember( MEMBER3, HOFFSET(h5_particle_struct, h5_particle_position_x), PredType::NATIVE_DOUBLE);
        mtype.insertMember( MEMBER4, HOFFSET(h5_particle_struct, h5_particle_position_y), PredType::NATIVE_DOUBLE);
        mtype.insertMember( MEMBER5, HOFFSET(h5_particle_struct, h5_particle_position_z), PredType::NATIVE_DOUBLE);

        // hsize_t dim[] = {5, particles_.size()};
        //std::cout << np << std::endl;
        hsize_t dim[] = {particles_.size()};
        DataSpace space(1, dim);

        //DataSet* dataset;
        //dataset = new DataSet(file_->createDataSet(DATASET_NAME, mtype, space));

        // create path.
        std::ostringstream ost_hdf5path;
        boost::scoped_ptr<Group> parent_group (new Group(file_->openGroup("/ParticleSpace")));
        ost_hdf5path << "/ParticleSpace/" << t;
        boost::scoped_ptr<Group> group (new Group(parent_group->createGroup( ost_hdf5path.str() )));

        std::string species_table_path = ost_hdf5path.str() + "/species";
        //std::string species_num_path = ost_hdf5path.str() + "/num";
        boost::scoped_ptr<H5::DataSet> dataset_id_table( new DataSet(file_->createDataSet(species_table_path, mtype, space)) );

        //dataset->write(h5_p, mtype);
        //dataset_id_table->write(particle_id_table.get(), mtype);
        dataset_id_table->write(h5_p, mtype);

    }

    std::pair<ParticleID, Particle> get_particle(ParticleID const& pid) const;
    std::vector<std::pair<ParticleID, Particle> > list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
    list_particles(Species const& species) const;

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius,
        ParticleID const& ignore) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        Position3 const& pos, Real const& radius,
        ParticleID const& ignore1, ParticleID const& ignore2) const;

private:

    void set_edge_lengths(Position3 const& edge_lengths);

protected:

    typedef utils::get_mapper_mf<
        ParticleID, particle_container_type::size_type>::type particle_map_type;

    Position3 edge_lengths_;
    particle_container_type particles_;
    particle_map_type index_map_;
};

} // ecell4

#endif /* __PARTICLE_SPACE_HPP */
