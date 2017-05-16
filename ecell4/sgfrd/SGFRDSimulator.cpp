#include "SGFRDSimulator.hpp"

namespace ecell4
{
namespace sgfrd
{
const Real SGFRDSimulator::single_circular_shell_factor        = 1.5;
const Real SGFRDSimulator::single_circular_shell_mergin        = 0.99;
const Real SGFRDSimulator::single_conical_surface_shell_factor = 1.5;
const Real SGFRDSimulator::single_conical_surface_shell_mergin = 0.99;

void SGFRDSimulator::domain_firer::operator()(const Single& dom)
{
    switch(dom.eventkind())
    {
        case Single::ESCAPE:
        {
            DUMP_MESSAGE("single escape");
            single_escapement executor(sim, dom);
            sim.execute_event(executor, dom.shell_id());
            return;
        }
        case Single::REACTION:
        {
            DUMP_MESSAGE("single reaction");
            single_reactor executor(sim, dom);
            sim.execute_event(executor, dom.shell_id());
            return;
        }
        case Single::UNKNOWN:
            throw std::logic_error("when firing Single: event unspecified");
        default:
            throw std::logic_error("when firing Single: invalid enum value");
    }
}

void SGFRDSimulator::domain_firer::operator()(const Pair& dom)
{//TODO
    std::cerr << "[WARNING] Pair domain firer has not been implemented yet"
              << std::endl;
    return;
}

void SGFRDSimulator::domain_firer::operator()(const Multi& dom)
{//TODO
    std::cerr << "[WARNING] Multi domain firer has not been implemented yet"
              << std::endl;
    return;
}

void SGFRDSimulator::single_burster::operator()(const circular_shell_type& sh)
{//TODO
    std::cerr << "bursting circualr shell" << std::endl;
    return;
}

void SGFRDSimulator::single_burster::operator()(const conical_surface_shell_type& sh)
{//TODO
    std::cerr << "bursting conical surface shell" << std::endl;
    return;
}

void SGFRDSimulator::pair_burster::operator()(const circular_shell_type& sh)
{//TODO
    std::cerr << "bursting pair circualr shell" << std::endl;
    return;
}

void SGFRDSimulator::pair_burster::operator()(const conical_surface_shell_type& sh)
{
    throw std::logic_error("bursting pair conical surface shell");
}


void SGFRDSimulator::single_escapement::operator()(const circular_shell_type& sh)
{
    DUMP_MESSAGE("single shell escapement circular shell");
    if(sh.size() == dom.particle().radius())
    {
        DUMP_MESSAGE("minimum shell. didnot move.");
        remnants.push_back(boost::make_tuple(dom.particle_id(), dom.particle(),
                                            sim.get_face_id(dom.particle_id())));
        return;
    }

    Particle   p   = dom.particle();
    ParticleID pid = dom.particle_id();

    const Real r   = sh.size() - p.radius();
    const Real theta = sim.uniform_real() * 2.0 * 3.141592653589793;
    DUMP_MESSAGE("r = " << r << ", theta = " << theta);
    const face_id_type   fid  = sim.get_face_id(pid);
    const triangle_type& face = sim.polygon().triangle_at(fid);
    const Real3 direction = rotate(theta, face.normal(), face.represent());
    DUMP_MESSAGE("direction = " << direction << ", length = " << length(direction));

    std::pair<std::pair<Real3, face_id_type>, Real3> state =
        std::make_pair(/*position = */std::make_pair(p.position(), fid),
                   /*displacement = */direction * r / length(direction));

    DUMP_MESSAGE("pos  = " << state.first.first << ", fid = " << state.first.second);
    unsigned int continue_count = 2;
    while(continue_count > 0)
    {
        state = sim.polygon().move_next_face(state.first, state.second);
        const Real3& disp = state.second;
        if(disp[0] == 0. && disp[1] == 0. && disp[2] == 0.) break;
        --continue_count;
        DUMP_MESSAGE("pos  = " << state.first.first << ", fid = " << state.first.second);
        DUMP_MESSAGE("disp = " << disp << ", length = " << length(disp));
    }
    if(continue_count == 0)
        std::cerr << "[WARNING] moving on face: precision lost" << std::endl;

    DUMP_MESSAGE("pos  = " << state.first.first << ", fid = " << state.first.second);
    DUMP_MESSAGE("disp = " << state.second << ", length = " << state.second);

    DUMP_MESSAGE("escaped.");

    p.position() = state.first.first;
    sim.update_particle(pid, p, state.first.second);
    remnants.push_back(boost::make_tuple(pid, p, state.first.second));
    return ;
}

void SGFRDSimulator::single_escapement::operator()(
        const conical_surface_shell_type& sh)
{
    Particle           p   = dom.particle();
    const ParticleID   pid = dom.particle_id();
    const face_id_type fid = sim.get_face_id(pid);
    DUMP_MESSAGE("escape-conical: pos  = " << p.position() << ", fid = " << fid);

    const Real r     = sh.size() - p.radius();
    greens_functions::GreensFunction2DRefWedgeAbs
        gf(p.D(), length(p.position() - sh.position()),
           r,     sh.shape().apex_angle());
    const Real theta = gf.drawTheta(sim.uniform_real(), r, dom.dt());

    DUMP_MESSAGE("escape-conical: r = " << r << ", theta = " << theta);

    const std::pair<Real3, face_id_type> state =
        sim.polygon().rotate_around_vertex(std::make_pair(p.position(), fid),
                                           sh.structure_id(), r, theta);

    DUMP_MESSAGE("escaped : pos = " << state.first << ", fid = " << state.second);

    p.position() = state.first;
    sim.update_particle(pid, p, state.second);
    remnants.push_back(boost::make_tuple(pid, p, state.second));
    return;
}

void SGFRDSimulator::single_reactor::operator()(const circular_shell_type& sh)
{//TODO
    std::cerr << "react in single circualr shell" << std::endl;
    return;
}

void SGFRDSimulator::single_reactor::operator()(const conical_surface_shell_type& sh)
{//TODO
    std::cerr << "react in single conical surface shell" << std::endl;
    return;
}

void SGFRDSimulator::create_event(
            const ParticleID& pid, const Particle& p, const face_id_type fid)
{
    DUMP_MESSAGE("create event");
    const std::pair<Real3, face_id_type> pos = std::make_pair(p.position(), fid);

    const Real min_circle_size = p.radius() * single_circular_shell_factor;
          Real max_circle_size = get_max_circle_size(pos);

    DUMP_MESSAGE("min circle size = " << min_circle_size);
    DUMP_MESSAGE("max circle size = " << max_circle_size);

    if(max_circle_size < min_circle_size)
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

        // TODO burst and form pair or multi if needed
        std::cerr << "[WARNING] intrusive domains exist." << std::endl;
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

Real SGFRDSimulator::distance_sq_to_segment(
        const Real3& p, const std::pair<Real3, Real3>& seg) const
{
    const Real3 ab = seg.second - seg.first;
    const Real3 ac = p - seg.first;
    const Real3 bc = p - seg.second;
    const Real dot = dot_product(ac, ab);
    if(dot <= 0.0) return length_sq(ac);
    const Real len = length_sq(ab);
    if(dot >= len) return length_sq(bc);
    return length_sq(ac) - (dot * dot) / len;
}

} // sgfrd
} // ecell4
