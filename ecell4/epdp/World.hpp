#ifndef WORLD_HPP
#define WORLD_HPP


#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/types.hpp>
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/Position3.hpp>
#include <ecell4/core/Context.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/core/extras.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/ParticleSpaceHDF5Writer.hpp>
#include <ecell4/core/Sphere.hpp>
#include "./ParticleTraits.hpp" // This refers ecell4::Particle
#include "structures.hpp"

#include "ParticleContainerBase.hpp"


#include <map>
#include <boost/lexical_cast.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/array.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/foreach.hpp>
#include "exceptions.hpp"
#include "generator.hpp"
#include "filters.hpp"
//#include "ParticleID.hpp"
//#include "SpeciesTypeID.hpp"
//#include "SpeciesInfo.hpp"
#include "Defs.hpp"
//#include "SerialIDGenerator.hpp"
#include "Transaction.hpp"
// #include "Structure.hpp"
// #include "Surface.hpp"
// #include "Region.hpp"
#include "geometry.hpp"
//#include "GSLRandomNumberGenerator.hpp"
//#include "Point.hpp" // XXX: workaround. should be removed later.
#include "Position3Type.hpp"
#include "utils/pair.hpp"

// #include "ParticleSimulationStructure.hpp"
// #include "CuboidalRegion.hpp"
// #include "PlanarSurface.hpp"
// #include "CylindricalSurface.hpp"
// #include "SphericalSurface.hpp"


// For twofold_container
inline
bool is_initialized(std::string const &obj)
{
    return (0 < obj.size());
}

template<typename Tderived_, typename TD_>
struct WorldTraitsBase
{
    typedef std::size_t size_type;
    typedef ecell4::Real length_type;
    typedef ecell4::Real D_type;
    typedef ecell4::Real time_type;
    typedef ecell4::ParticleID particle_id_type;
    typedef ecell4::SerialIDGenerator<particle_id_type> particle_id_generator;
    typedef ecell4::Species::serial_type species_id_type; // std::string
    typedef ecell4::Particle particle_type;
    typedef ecell4::Position3 position_type;
    typedef ecell4::GSLRandomNumberGenerator rng_type;
    typedef ecell4::Model model_type;

    struct MoleculeInfo
    {
        const ecell4::Real radius;
        const ecell4::Real D;
        const std::string structure_id;
    };

    typedef MoleculeInfo molecule_info_type;
    // typedef MoleculeInfo species_info_type;
    // typedef SpeciesInfo<species_id_type, D_type, length_type, structure_id_type>
    //     species_info_type;

    // typedef Sphere particle_shape_type;
    typedef ecell4::Sphere particle_shape_type;
    typedef std::string structure_id_type;

    typedef std::pair<particle_id_type, particle_type> particle_id_pair;
    // typedef std::pair<const particle_id_type, particle_type> particle_id_pair;
    typedef std::pair<particle_id_pair, length_type> particle_id_pair_and_distance;
    // typedef unassignable_adapter<particle_id_pair_and_distance, get_default_impl::std::vector> particle_id_pair_and_distance_list;
    typedef std::vector<particle_id_pair_and_distance> particle_id_pair_and_distance_list;
    typedef abstract_limited_generator<particle_id_pair> particle_id_pair_generator;

    typedef ecell4::Structure<Tderived_> structure_type;
    typedef ecell4::Structure<Tderived_> particle_simulation_structure_type;
    typedef ecell4::AABBRegion<Tderived_> cuboidal_region_type;

    // typedef Structure<Tderived_> structure_type;
    // typedef ParticleSimulationStructure<Tderived_>
    //     particle_simulation_structure_type;
    // // typedef Surface<Tderived_> surface_type;
    // // typedef Region<Tderived_> region_type;
    // // typedef SphericalSurface<Tderived_> spherical_surface_type;
    // // typedef CylindricalSurface<Tderived_> cylindrical_surface_type;
    // // typedef PlanarSurface<Tderived_> planar_surface_type;
    // typedef CuboidalRegion<Tderived_> cuboidal_region_type;

    static const Real TOLERANCE = 1e-7;
};

template<typename TD_>
struct WorldTraits: public WorldTraitsBase<WorldTraits<TD_>, TD_>
{
public:
    typedef WorldTraitsBase<WorldTraits<TD_>, TD_> base_type;
    typedef typename base_type::length_type length_type;
    typedef typename base_type::position_type position_type;

    template<typename Tval_>
    static Tval_ apply_boundary(Tval_ const& v, position_type const& edge_lengths)
    {
        return v;
    }

    template<typename Tval_>
    static Tval_ cyclic_transpose(Tval_ const& p0, Tval_ const& p1, Tval_ const& world_size)
    {
        return p0;
    }

    template<typename T1_, typename T2_>
    static length_type distance(T1_ const& p0, T2_ const& p1, position_type const& edge_lengths)
    {
        return ::distance(p0, p1);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void each_neighbor(Toc_& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor(oc.index(pos), fun);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void each_neighbor(Toc_ const& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor(oc.index(pos), fun);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor(oc, fun, cmp);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_ const& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor(oc, fun, cmp);
    }
};

template<typename TD_>
struct CyclicWorldTraits: public WorldTraitsBase<CyclicWorldTraits<TD_>, TD_>
{
public:
    typedef WorldTraitsBase<CyclicWorldTraits<TD_>, TD_> base_type;
    typedef typename base_type::length_type length_type;
    typedef typename base_type::position_type position_type;

    template<typename Tval_>
    static Tval_ apply_boundary(Tval_ const& v, position_type const& edge_lengths)
    {
        return ::apply_boundary(v, edge_lengths);
    }

    static length_type cyclic_transpose(length_type const& p0, length_type const& p1, length_type const& world_size)
    {
        return ::cyclic_transpose(p0, p1, world_size);
    }

    static position_type cyclic_transpose(position_type const& p0, position_type const& p1, position_type const& edge_lengths)
    {
        return ::cyclic_transpose(p0, p1, edge_lengths);
    }

    template<typename T1_, typename T2_, typename T3_>
    static length_type distance(T1_ const& p0, T2_ const& p1, T3_ const& edge_lengths)
    {
        return distance_cyclic(p0, p1, edge_lengths);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void each_neighbor(Toc_& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor_cyclic(oc.index(pos), fun);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void each_neighbor(Toc_ const& oc, Tfun_& fun, Tsphere_ const& pos)
    {
        oc.each_neighbor_cyclic(oc.index(pos), fun);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor_cyclic(oc, fun, cmp);
    }

    template<typename Toc_, typename Tfun_, typename Tsphere_>
    static void take_neighbor(Toc_ const& oc, Tfun_& fun, const Tsphere_& cmp)
    {
        take_neighbor_cyclic(oc, fun, cmp);
    }
};

template<typename Ttraits_>
class World: public ParticleContainerBase<World<Ttraits_>, Ttraits_>
{
public:
    typedef Ttraits_ traits_type;
    typedef ParticleContainerBase<World> base_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef typename traits_type::length_type length_type;
    typedef typename traits_type::molecule_info_type molecule_info_type;
    typedef typename traits_type::position_type position_type;
    typedef typename traits_type::particle_type particle_type;
    typedef typename traits_type::particle_id_type particle_id_type;
    typedef typename traits_type::particle_id_generator particle_id_generator;
    typedef typename traits_type::species_id_type species_id_type;
    typedef typename traits_type::particle_shape_type particle_shape_type;
    typedef typename traits_type::size_type size_type;
    typedef typename traits_type::structure_id_type structure_id_type;
    typedef typename traits_type::structure_type structure_type;
    typedef typename traits_type::rng_type rng_type;
    typedef typename traits_type::particle_id_pair particle_id_pair;
    typedef typename traits_type::particle_id_pair_and_distance_list
        particle_id_pair_and_distance_list;
    typedef typename traits_type::model_type model_type;

    typedef typename base_type::matrix_sizes_type matrix_sizes_type;

protected:
    typedef std::map<species_id_type, molecule_info_type> molecule_info_map;
    typedef std::map<structure_id_type, boost::shared_ptr<structure_type> > structure_map;
    typedef std::set<particle_id_type> particle_id_set;
    typedef std::map<species_id_type, particle_id_set> per_species_particle_id_set;
    typedef select_second<typename molecule_info_map::value_type> species_second_selector_type;
    typedef select_second<typename structure_map::value_type> surface_second_selector_type;

public:
    typedef boost::transform_iterator<species_second_selector_type,
            typename molecule_info_map::const_iterator> molecule_info_iterator;
    typedef boost::transform_iterator<surface_second_selector_type,
            typename structure_map::const_iterator> surface_iterator;
    typedef sized_iterator_range<molecule_info_iterator> molecule_info_range;
    typedef sized_iterator_range<surface_iterator> structures_range;

public:

    World(
        const position_type& edge_lengths = position_type(1, 1, 1),
        const matrix_sizes_type& sizes = matrix_sizes_type(3, 3, 3))
        : base_type(edge_lengths, sizes)
    {
        rng_ = boost::shared_ptr<rng_type>(new rng_type());
        (*rng_).seed();

        add_world_structure();
    }

    World(
        const position_type& edge_lengths, const matrix_sizes_type& sizes,
        const boost::shared_ptr<rng_type>& rng)
        : base_type(edge_lengths, sizes), rng_(rng)
    {
        add_world_structure();
    }

    World(
        const position_type& edge_lengths, const matrix_sizes_type& sizes,
        const boost::shared_ptr<ecell4::RandomNumberGenerator>& rng)
        : base_type(edge_lengths, sizes),
        rng_(boost::dynamic_pointer_cast<rng_type>(rng))
    {
        add_world_structure();
    }

    World(const std::string filename)
        : base_type(position_type(1, 1, 1), matrix_sizes_type(3, 3, 3)), rng_()
    {
        rng_ = boost::shared_ptr<rng_type>(new rng_type());
        this->load(filename);
    }

    virtual particle_id_pair new_particle(species_id_type const& sid,
            position_type const& pos)
    {
        molecule_info_type const& minfo(get_molecule_info(sid));
        particle_id_pair retval(pidgen_(),
            particle_type(sid, pos, minfo.radius, minfo.D));
        update_particle(retval);
        return retval;
    }

    virtual bool update_particle(particle_id_pair const& pi_pair)
    {
        std::pair<bool, typename base_type::particle_matrix_type::iterator>
            retval(base_type::__has_particle(pi_pair.first));
        if (retval.first)
        {
            typename base_type::particle_matrix_type::iterator&
                i(retval.second);
            if ((*i).second.sid() != pi_pair.second.sid())
            {
                particle_pool_[(*i).second.sid()].erase((*i).first);

                typename per_species_particle_id_set::iterator
                    j(particle_pool_.find(pi_pair.second.sid()));
                if (j == particle_pool_.end())
                {
                    this->register_species(pi_pair);
                    // this->register_species(pi_pair.second.sid());
                    j = particle_pool_.find(pi_pair.second.sid());
                }
                (*j).second.insert(pi_pair.first);
            }
            base_type::__update_particle(i, pi_pair);
            return false;
        }

        const bool is_succeeded(base_type::update_particle(pi_pair));
        BOOST_ASSERT(is_succeeded);
        // BOOST_ASSERT(base_type::update_particle(pi_pair)); //XXX: -DNDEBUG!!!
        typename per_species_particle_id_set::iterator
            k(particle_pool_.find(pi_pair.second.sid()));
        if (k == particle_pool_.end())
        {
            this->register_species(pi_pair);
            // this->register_species(pi_pair.second.sid());
            k = particle_pool_.find(pi_pair.second.sid());
        }
        (*k).second.insert(pi_pair.first);
        return true;
    }

    virtual bool remove_particle(particle_id_type const& id)
    {
        bool found(false);
        particle_id_pair pp(this->get_particle(id, found));
        if (!found)
        {
            return false;
        }
        particle_pool_[pp.second.sid()].erase(id);
        base_type::remove_particle(id);
        return true;
    }

    virtual molecule_info_type const& get_molecule_info(species_id_type const& sid)
    {
        typename molecule_info_map::const_iterator i(molecule_info_map_.find(sid));
        if (molecule_info_map_.end() == i)
        {
            return this->register_species(sid);
            // throw not_found(std::string("Unknown species (id=")
            //     + boost::lexical_cast<std::string>(sid) + ")");
        }
        return (*i).second;
    }

    virtual molecule_info_type const& find_molecule_info(species_id_type const& sid) const
    {
        typename molecule_info_map::const_iterator i(molecule_info_map_.find(sid));
        if (molecule_info_map_.end() == i)
        {
            throw not_found(std::string("Unknown species (id=")
                + boost::lexical_cast<std::string>(sid) + ")");
        }
        return (*i).second;
    }

    molecule_info_range get_molecule_info_range() const
    {
        return molecule_info_range(
            molecule_info_iterator(
                molecule_info_map_.begin(), species_second_selector_type()),
            molecule_info_iterator(
                molecule_info_map_.end(), species_second_selector_type()),
            molecule_info_map_.size());
    }

    bool add_structure(boost::shared_ptr<structure_type> surface)
    {
        return structure_map_.insert(std::make_pair(surface->id(), surface)).second;
    }

    virtual boost::shared_ptr<structure_type> get_structure(
        structure_id_type const& id) const
    {
        typename structure_map::const_iterator i(structure_map_.find(id));
        if (structure_map_.end() == i)
        {
            throw not_found(std::string("Unknown surface (id=")
                + boost::lexical_cast<std::string>(id) + ")");
        }
        return (*i).second;
    }

    structures_range get_structures() const
    {
        return structures_range(
            surface_iterator(structure_map_.begin(), surface_second_selector_type()),
            surface_iterator(structure_map_.end(), surface_second_selector_type()),
            structure_map_.size());
    }

    particle_id_set get_particle_ids(species_id_type const& sid) const
    {
        typename per_species_particle_id_set::const_iterator i(
            particle_pool_.find(sid));
        if (i == particle_pool_.end())
        {
            throw not_found(std::string("Unknown species (id=")
                + boost::lexical_cast<std::string>(sid) + ")");
        }
        return (*i).second;
    }

    /** ecell4::Space
     */

    inline boost::shared_ptr<rng_type>& rng()
    {
        return rng_;
    }

    virtual void save(const std::string& filename) const
    {
        boost::scoped_ptr<H5::H5File>
            fout(new H5::H5File(filename.c_str(), H5F_ACC_TRUNC));
        rng_->save(fout.get());
        pidgen_.save(fout.get());
        boost::scoped_ptr<H5::Group>
            group(new H5::Group(fout->createGroup("ParticleSpace")));
        //  ps_->save(group.get());
        ecell4::save_particle_space(*this, group.get());

        /** matrix_sizes
         */
        const matrix_sizes_type sizes = base_type::matrix_sizes();
        const hsize_t dims[] = {3};
        const H5::ArrayType sizes_type(H5::PredType::NATIVE_INT, 1, dims);
        H5::Attribute attr_sizes(
            group->createAttribute(
                "matrix_sizes", sizes_type, H5::DataSpace(H5S_SCALAR)));
        int data[] = {sizes[0], sizes[1], sizes[2]};
        attr_sizes.write(sizes_type, data);
    }

    virtual void load(const std::string& filename)
    {
        //XXX: structures will be lost.
        //XXX: the order of particles in MatrixSpace will be lost.
        //XXX: initialize Simulator
        boost::scoped_ptr<H5::H5File>
            fin(new H5::H5File(filename.c_str(), H5F_ACC_RDONLY));
        const H5::Group group(fin->openGroup("ParticleSpace"));

        /** matrix_sizes
         */
        int data[3];
        const hsize_t dims[] = {3};
        const H5::ArrayType sizes_type(H5::PredType::NATIVE_INT, 1, dims);
        group.openAttribute("matrix_sizes").read(sizes_type, data);
        matrix_sizes_type sizes(data[0], data[1], data[2]);
        //XXX: reset is called twice. see ecell4::load_particle_space
        this->reset(edge_lengths(), sizes);

        // ps_->load(group);
        ecell4::load_particle_space(group, this);
        pidgen_.load(*fin);
        rng_->load(*fin);
    }

    virtual void clear()
    {
        // particle_id_generator pidgen_;
        // boost::shared_ptr<rng_type> rng_;
        // boost::weak_ptr<model_type> model_;
        ; // do nothing

        // molecule_info_map molecule_info_map_;
        // structure_map structure_map_;
        // per_species_particle_id_set particle_pool_;
        molecule_info_map_.clear();
        structure_map_.clear();
        particle_pool_.clear();

        // ParticleContainerBase
        // particle_matrix_type pmat_;
        // time_type t_;
        base_type::clear();
    }

    virtual const length_type volume() const
    {
        const position_type& L(edge_lengths());
        return L[0] * L[1] * L[2];
    }

    virtual const position_type& edge_lengths() const
    {
        return base_type::edge_lengths();
    }

    virtual void reset(const position_type& lengths)
    {
        base_type::reset(lengths);
    }

    virtual void reset(const position_type& lengths, const matrix_sizes_type& sizes)
    {
        base_type::reset(lengths, sizes);
    }

    virtual bool has_particle(const particle_id_type& pid) const
    {
        return base_type::has_particle(pid);
    }

    virtual ecell4::Integer num_particles() const
    {
        return base_type::num_particles();
    }

    virtual ecell4::Integer num_particles_exact(const ecell4::Species& sp) const
    {
        typename per_species_particle_id_set::const_iterator
            i(particle_pool_.find(sp.serial()));
        if (i == particle_pool_.end())
        {
            return 0;
        }
        return (*i).second.size();
    }

    virtual ecell4::Integer num_particles(const ecell4::Species& sp) const
    {
        ecell4::Integer retval(0);
        ecell4::SpeciesExpressionMatcher sexp(sp);
        for (typename per_species_particle_id_set::const_iterator
            i(particle_pool_.begin()); i != particle_pool_.end(); ++i)
        {
            const ecell4::Species tgt((*i).first);
            if (sexp.match(tgt))
            {
                retval += (*i).second.size();
            }
        }
        return retval;
    }

    virtual ecell4::Real get_value(const ecell4::Species& sp) const
    {
        return static_cast<ecell4::Real>(num_molecules(sp));
    }

    virtual ecell4::Real get_value_exact(const ecell4::Species& sp) const
    {
        return static_cast<ecell4::Real>(num_molecules_exact(sp));
    }

    virtual ecell4::Integer num_molecules(const ecell4::Species& sp) const
    {
        ecell4::Integer retval(0);
        ecell4::SpeciesExpressionMatcher sexp(sp);
        for (typename per_species_particle_id_set::const_iterator
            i(particle_pool_.begin()); i != particle_pool_.end(); ++i)
        {
            const ecell4::Species tgt((*i).first);
            retval += sexp.count(tgt) * (*i).second.size();
        }
        return retval;
    }

    virtual ecell4::Integer num_molecules_exact(const ecell4::Species& sp) const
    {
        return num_particles_exact(sp);
    }

    virtual ecell4::Integer num_species() const
    {
        return particle_pool_.size();
    }

    virtual bool has_species(const ecell4::Species& sp) const
    {
        return (particle_pool_.find(sp.serial()) != particle_pool_.end());
    }

    virtual std::vector<std::pair<particle_id_type, particle_type> > list_particles() const
    {
        std::vector<std::pair<particle_id_type, particle_type> > retval;
        retval.reserve(num_particles());
        BOOST_FOREACH(particle_id_pair p, this->get_particles_range())
        {
            retval.push_back(p);
        }
        return retval;
    }

    virtual std::vector<std::pair<particle_id_type, particle_type> >
        list_particles(const ecell4::Species& sp) const
    {
        std::vector<std::pair<particle_id_type, particle_type> > retval;
        ecell4::SpeciesExpressionMatcher sexp(sp);
        for (typename per_species_particle_id_set::const_iterator
            i(particle_pool_.begin()); i != particle_pool_.end(); ++i)
        {
            const ecell4::Species tgt((*i).first);
            if (sexp.match(tgt))
            {
                for (typename particle_id_set::const_iterator j((*i).second.begin());
                    j != (*i).second.end(); ++j)
                {
                    const particle_id_type& pid(*j);
                    retval.push_back(this->get_particle(pid));
                }
            }
        }
        return retval;
    }

    virtual std::vector<std::pair<particle_id_type, particle_type> >
        list_particles_exact(const ecell4::Species& sp) const
    {
        std::vector<std::pair<particle_id_type, particle_type> > retval;
        typename per_species_particle_id_set::const_iterator
            i(particle_pool_.find(sp.serial()));
        if (i == particle_pool_.end())
        {
            return retval;
        }

        for (typename particle_id_set::const_iterator j((*i).second.begin());
            j != (*i).second.end(); ++j)
        {
            const particle_id_type& pid(*j);
            retval.push_back(this->get_particle(pid));
        }
        return retval;
    }

    std::vector<std::pair<std::pair<particle_id_type, particle_type>, length_type> >
    list_particles_within_radius(
        const position_type& pos, const length_type& radius) const
    {
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
            base_type::check_overlap(particle_shape_type(pos, radius)));
        std::vector<std::pair<std::pair<particle_id_type, particle_type>, length_type> >
            retval;
        if (overlapped && ::size(*overlapped))
        {
            for (typename particle_id_pair_and_distance_list::const_iterator
                i(overlapped->begin()); i != overlapped->end(); ++i)
            {
                retval.push_back(*i);
            }
        }
        return retval;
    }

    std::vector<std::pair<std::pair<particle_id_type, particle_type>, length_type> >
    list_particles_within_radius(
        const position_type& pos, const length_type& radius,
        const particle_id_type& ignore) const
    {
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
            base_type::check_overlap(particle_shape_type(pos, radius), ignore));
        std::vector<std::pair<std::pair<particle_id_type, particle_type>, length_type> >
            retval;
        if (overlapped && ::size(*overlapped))
        {
            for (typename particle_id_pair_and_distance_list::const_iterator
                i(overlapped->begin()); i != overlapped->end(); ++i)
            {
                retval.push_back(*i);
            }
        }
        return retval;
    }

    std::vector<std::pair<std::pair<particle_id_type, particle_type>, length_type> >
    list_particles_within_radius(
        const position_type& pos, const length_type& radius,
        const particle_id_type& ignore1, const particle_id_type& ignore2) const
    {
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
            base_type::check_overlap(
                particle_shape_type(pos, radius), ignore1, ignore2));
        std::vector<std::pair<std::pair<particle_id_type, particle_type>, length_type> >
            retval;
        if (overlapped && ::size(*overlapped))
        {
            for (typename particle_id_pair_and_distance_list::const_iterator
                i(overlapped->begin()); i != overlapped->end(); ++i)
            {
                retval.push_back(*i);
            }
        }
        return retval;
    }

    void bind_to(boost::shared_ptr<model_type> model)
    {
        if (boost::shared_ptr<model_type> bound_model = lock_model())
        {
            if (bound_model.get() != model.get())
            {
                std::cerr << "Warning: Model already bound to BDWorld"
                    << std::endl;
            }
        }
        model_ = model;
    }

    boost::shared_ptr<model_type> lock_model() const
    {
        return model_.lock();
    }

    /**
     * This is a function in the traits of ecell4::ParticleSpace.
     * Be carefull about the difference from
     * "particle_id_pair new_particle(species_id_type const&, position_type const&)".
     */
    std::pair<std::pair<particle_id_type, particle_type>, bool>
    new_particle(const ecell4::Species& sp, const position_type& pos)
    {
        const species_id_type sid(sp.serial());
        molecule_info_type const& minfo(get_molecule_info(sid));
        return new_particle(particle_type(sid, pos, minfo.radius, minfo.D));
    }

    std::pair<std::pair<particle_id_type, particle_type>, bool>
    new_particle(const particle_type& p)
    {
        particle_id_pair retval(pidgen_(), p);
        boost::scoped_ptr<particle_id_pair_and_distance_list> overlapped(
            base_type::check_overlap(
                particle_shape_type(p.position(), p.radius())));
        if (overlapped && ::size(*overlapped))
        {
            return std::make_pair(retval, false);
        }
        else
        {
            return std::make_pair(retval, update_particle(retval));
        }
    }

    virtual bool update_particle(const particle_id_type& pid, const particle_type& p)
    {
        return this->update_particle(std::make_pair(pid, p));
    }

    void add_molecules(const ecell4::Species& sp, const ecell4::Integer& num)
    {
        ecell4::extras::throw_in_particles(*this, sp, num, rng());
    }

    void add_molecules(
        const ecell4::Species& sp, const ecell4::Integer& num,
        const ecell4::Shape& shape)
    {
        ecell4::extras::throw_in_particles(*this, sp, num, shape, rng());
    }

    void remove_molecules(const ecell4::Species& sp, const ecell4::Integer& num)
    {
        if (num == 0)
        {
            return;
        }
        else if (num < 0)
        {
            throw std::invalid_argument(
                "The number of molecules must be positive.");
        }

        typename per_species_particle_id_set::const_iterator
            i(particle_pool_.find(sp.serial()));
        if (i == particle_pool_.end() || (*i).second.size() < num)
        {
            throw std::invalid_argument(
                "The number of molecules cannot be negative.");
        }

        for (unsigned int j(0); j < num; ++j)
        {
            const Integer n(rng()->uniform_int(0, (*i).second.size() - 1));
            // typename particle_id_set::const_iterator
            //     target(std::advance((*i).second.begin(), n));
            // this->remove_particle((*target).first);
            typename particle_id_set::const_iterator
                target((*i).second.begin());
            std::advance(target, n);
            this->remove_particle(*target);
        }
    }

    /**
     * draw attributes of species and return it as a molecule info.
     * @param sp a species
     * @return info a molecule info
     */
    molecule_info_type get_molecule_info(const ecell4::Species& sp) const
    {
        const molecule_info_type defaults = {0.0, 0.0, "world"};
        return get_molecule_info(sp, defaults);
    }

    molecule_info_type get_molecule_info(
        const ecell4::Species& sp, const molecule_info_type& defaults) const
    {
        ecell4::Real radius(defaults.radius), D(defaults.D);
        std::string structure_id(defaults.structure_id);

        if (sp.has_attribute("radius") && sp.has_attribute("D"))
        {
            radius = std::atof(sp.get_attribute("radius").c_str());
            D = std::atof(sp.get_attribute("D").c_str());
            if (sp.has_attribute("structure_id"))
            {
                structure_id = sp.get_attribute("structure_id");
            }
        }
        else if (boost::shared_ptr<model_type> bound_model = lock_model())
        {
            ecell4::Species attributed(bound_model->apply_species_attributes(sp));

            if (attributed.has_attribute("radius")
                && attributed.has_attribute("D"))
            {
                radius = std::atof(
                    attributed.get_attribute("radius").c_str());
                D = std::atof(attributed.get_attribute("D").c_str());
            }

            if (sp.has_attribute("structure_id"))
            {
                structure_id = attributed.get_attribute("structure_id");
            }
        }

        molecule_info_type info = {radius, D, structure_id};
        return info;
    }

protected:

    const molecule_info_type& register_species(const particle_id_pair& pid_pair)
    {
        const molecule_info_type defaults
            = {pid_pair.second.radius(), pid_pair.second.D(), "world"};
        const species_id_type sid(pid_pair.second.sid());
        const ecell4::Species sp(sid);
        molecule_info_type info(get_molecule_info(sp, defaults));
        molecule_info_map_.insert(std::make_pair(sid, info));
        particle_pool_[sid] = particle_id_set();
        return (*molecule_info_map_.find(sid)).second;
    }

    const molecule_info_type& register_species(const species_id_type& sid)
    {
        const ecell4::Species sp(sid);
        molecule_info_type info(get_molecule_info(sp));
        molecule_info_map_.insert(std::make_pair(sid, info));
        particle_pool_[sid] = particle_id_set();
        return (*molecule_info_map_.find(sid)).second;
    }

    void add_world_structure()
    {
        typedef typename traits_type::cuboidal_region_type cuboidal_region_type;
        typedef typename cuboidal_region_type::shape_type
            cuboidal_region_shape_type;

        this->add_structure(
            boost::shared_ptr<structure_type>(
                new cuboidal_region_type(
                    "world", cuboidal_region_shape_type(
                        position_type(0, 0, 0), edge_lengths()))));
        // const position_type& center(edge_lengths() * 0.5);
        // this->add_structure(
        //     boost::shared_ptr<structure_type>(
        //         new cuboidal_region_type(
        //             "world", cuboidal_region_shape_type(center, center))));
    }

private:
    particle_id_generator pidgen_;
    molecule_info_map molecule_info_map_;
    structure_map structure_map_;
    per_species_particle_id_set particle_pool_;

    /** ecell4::Space
     */
    boost::shared_ptr<rng_type> rng_;
    boost::weak_ptr<model_type> model_;
};

#endif /* WORLD_HPP */
