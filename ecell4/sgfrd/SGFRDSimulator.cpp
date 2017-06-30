#include "SGFRDSimulator.hpp"

namespace ecell4
{
namespace sgfrd
{
const Real SGFRDSimulator::single_circular_shell_factor        = 1.5;
const Real SGFRDSimulator::single_circular_shell_mergin        = 1.0 - 1e-7;
const Real SGFRDSimulator::single_conical_surface_shell_factor = 1.5;
const Real SGFRDSimulator::single_conical_surface_shell_mergin = 1.0 - 1e-7;

void SGFRDSimulator::fire_single(const Single& dom, DomainID did)
{
    const ShellID sid(dom.shell_id());
    ParticleID pid; Particle p; FaceID fid;
    switch(dom.eventkind())
    {
    case Single::ESCAPE:
    {
        DUMP_MESSAGE("single escape");
        boost::tie(pid, p, fid) = boost::apply_visitor(make_visitor(
            resolve<const circular_shell_type&,
                    boost::tuple<ParticleID, Particle, FaceID> >(boost::bind(
                &self_type::escape_single<circular_shell_type>,
                this, _1, dom)),
            resolve<const conical_surface_shell_type&,
                    boost::tuple<ParticleID, Particle, FaceID> >(boost::bind(
                &self_type::escape_single<conical_surface_shell_type>,
                this, _1, dom))
            ), get_shell(sid));
        this->remove_shell(sid);
        this->create_event(pid, p, fid);
        return;
    }
    case Single::REACTION:
    {
        DUMP_MESSAGE("single reaction");
        BOOST_AUTO(results, boost::apply_visitor(make_visitor(
            resolve<const circular_shell_type&,
                    boost::container::static_vector<
                        boost::tuple<ParticleID, Particle, FaceID>, 2>
                    >(boost::bind(
                &self_type::reaction_single<circular_shell_type>,
                this, _1, dom, did)),
            resolve<const conical_surface_shell_type&,
                    boost::container::static_vector<
                        boost::tuple<ParticleID, Particle, FaceID>, 2>
                    >(boost::bind(
                &self_type::reaction_single<conical_surface_shell_type>,
                this, _1, dom, did))
            ), get_shell(sid)));
        this->remove_shell(sid);

        BOOST_FOREACH(boost::tie(pid, p, fid), results)
        {
            this->create_event(pid, p, fid);
        }
        return;
    }
    case Single::UNKNOWN:
        throw std::logic_error("when firing Single: event unspecified");
    default:
        throw std::logic_error("when firing Single: invalid enum value");
    }
}


SGFRDSimulator::bursted_type
SGFRDSimulator::burst_overlaps(const Particle& p, const face_id_type& fid)
{
    const Real tm = this->time();
    BOOST_AUTO(intruders, this->get_intrusive_domains(
                std::make_pair(p.position(), fid), p.radius()));

    bursted_type bursted;
    DomainID did;
    BOOST_FOREACH(boost::tie(did, boost::tuples::ignore), intruders)
    {
        BOOST_FOREACH(typename bursted_type::value_type const& elem,
                      burst_event(std::make_pair(did, pickout_event(did)), tm))
        {
            bursted.push_back(elem);
        }
    }
    return bursted;
}

DomainID SGFRDSimulator::create_event(
            const ParticleID& pid, const Particle& p, const face_id_type fid)
{
    DUMP_MESSAGE("create event");
    const std::pair<Real3, face_id_type> pos = std::make_pair(p.position(), fid);

    const Real min_circle_size = p.radius() * single_circular_shell_factor;
          Real max_circle_size = get_max_circle_size(pos);

    DUMP_MESSAGE("min circle size = " << min_circle_size);
    DUMP_MESSAGE("max circle size = " << max_circle_size);

    if(max_circle_size < min_circle_size)// draw conical shell
    {
        DUMP_MESSAGE("drawing conical shell: min_circle_size = " << min_circle_size);
        const std::vector<std::pair<vertex_id_type, Real> > intrusive_vertices(
                get_intrusive_vertices(pos, std::numeric_limits<Real>::infinity()));

        const vertex_id_type& vid = intrusive_vertices.front().first;
        DUMP_MESSAGE("vertex id = " << vid << ", distance = "
                     << intrusive_vertices.front().second);

        const Real min_cone_size =
            (p.radius() + intrusive_vertices.front().second) *
            single_conical_surface_shell_factor;
        const Real max_cone_size = get_max_cone_size(vid);
        DUMP_MESSAGE("min cone size = " << min_cone_size);
        DUMP_MESSAGE("max cone size = " << max_cone_size);

        const std::vector<std::pair<DomainID, Real> > intrusive_domains(
                get_intrusive_domains(vid, max_cone_size));
        DUMP_MESSAGE("intrusive domains = " << intrusive_domains.size());

        if(intrusive_domains.empty())
        {
            return add_event(create_single(create_single_conical_surface_shell(
                                           vid, max_cone_size), pid, p));
        }

        if(intrusive_domains.front().second > min_cone_size)
        {
            DUMP_MESSAGE("avoid domains overlapping: distance = "
                         << intrusive_domains.front().second);
            return add_event(create_single(create_single_conical_surface_shell(
                vid, intrusive_domains.front().second *
                single_conical_surface_shell_mergin), pid, p));
        }

        //TODO: burst intruder_domains and get new positions of particles
        //note: L3347-3367 in egfrd/EGFRDSimulator.hpp
        DUMP_MESSAGE("[WARNING] ignoring intrusive domains.");
        return add_event(create_single(create_single_conical_surface_shell(
                                       vid, max_cone_size), pid, p));
    }

    DUMP_MESSAGE("drawing circular shell");

    const std::vector<std::pair<DomainID, Real> > intrusive_domains(
            get_intrusive_domains(pos, max_circle_size));
    DUMP_MESSAGE("intrusive domains: " << intrusive_domains.size());

    if(intrusive_domains.empty())
    {
        return add_event(create_single(create_single_circular_shell(
                         pos, max_circle_size), pid, p));
    }

    if(intrusive_domains.front().second > min_circle_size)
    {
        return add_event(create_single(create_single_circular_shell(
            pos, intrusive_domains.front().second * single_circular_shell_mergin),
            pid, p));
    }

    // TODO burst!
    std::cerr << "[WARNING] ignoring intrusive domains." << std::endl;
    return add_event(create_single(create_single_circular_shell(
                     pos, max_circle_size), pid, p));
}

} // sgfrd
} // ecell4
