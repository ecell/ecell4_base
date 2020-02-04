#include <ecell4/sgfrd/SGFRDWorld.hpp>
#include <boost/container/static_vector.hpp>

namespace ecell4
{
namespace sgfrd
{

std::pair<std::pair<ParticleID, Particle>, bool>
SGFRDWorld::new_particle(const Particle& p)
{
    if(const auto pfid = this->find_face(p.position()))
    {
        Particle p_(p);
        p_.position() = pfid->first;
        return this->new_particle(p_, pfid->second);
    }
    throw std::invalid_argument("[error] SGFRDWorld::new_particle: "
            "particle locates distant from polygon");
}

std::pair<std::pair<ParticleID, Particle>, bool>
SGFRDWorld::new_particle(const Particle& p, const FaceID& fid)
{
    if(p.radius() > this->estimated_possible_largest_particle_radius_)
    {
        throw NotSupported("[error] Particle size exceeds the estimated limit. "
                "particle size must be smaller than the width of triangles");
    }

    const ParticleID pid = pidgen_();
    // now this consider only 2D particles
    const auto overlap2d(list_particles_within_radius(
            std::make_pair(p.position(), fid), p.radius()));

    if(!overlap2d.empty())
    {
//         std::cout << "overlapped particle = {";
//         std::cout << '{' << overlap2d.front().first.first
//                   << ", distance = " << overlap2d.front().second << '}';
//         for(std::size_t i=1; i<overlap2d.size(); ++i)
//         {
//             std::cout << ", {" << overlap2d.at(i).first.first
//                       << ", distance = " << overlap2d.at(i).second << '}';
//         }
//         std::cout << '}' << std::endl;
        return std::make_pair(std::make_pair(pid, p), false);
    }
    return std::make_pair(std::make_pair(pid, p), update_particle(pid, p, fid));
}

std::pair<std::pair<ParticleID, Particle>, bool>
SGFRDWorld::throw_in_particle(const Species& sp)
{
    const auto info = this->get_molecule_info(sp);

    Real3 pos; FaceID fid;
    pos = this->polygon_->draw_position(this->rng_, fid);
    const Particle p(sp, pos, info.radius, info.D);
    return this->new_particle(p, fid);
}

void SGFRDWorld::add_molecules(const Species& sp, const Integer& num)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }
    for(Integer i=0; i<num; ++i)
    {
        while(this->throw_in_particle(sp).second == false)
        {
            /*do nothing*/
        }
    }
    return;
}

void SGFRDWorld::add_molecules(const Species& sp, const Integer& num,
                               const boost::shared_ptr<Shape> shape)
{
    if (num < 0)
    {
        throw std::invalid_argument("The number of molecules must be positive.");
    }
    else if (num == 0)
    {
        return;
    }

    // XXX: this implementation is not only inefficient, but also unsafe because
    //      it may not stop. Since `shape` is a base class, here the concrete
    //      representation of a shape cannot be obtained (except dynamic_cast).
    //      There is no way to predict the precise area of the overlapped
    //      region. If the overlapped region is just a point or a line without
    //      area, no particle can be inside it. If the overlapped region is too
    //      small to place `num` number of particles, it also does not finish.
    //      But still some faces overlap with the shape, and there is no way to
    //      make it sure that this function never ends, we need to continue
    //      searching.

    const auto& width = this->edge_lengths();
    const auto  transpose_direction = [&width](const Triangle& triangle)
        noexcept -> boost::container::static_vector<Real3, 3> {
            const auto& v1 = triangle.vertices()[0];
            const auto& v2 = triangle.vertices()[1];
            const auto& v3 = triangle.vertices()[2];

            boost::container::static_vector<Real3, 3> retval;
            for(std::size_t i=0; i<3; ++i)
            {
                Real3 disp(0.0, 0.0, 0.0);
                if(v1[i] < 0.0 || v2[i] < 0.0 || v3[i] < 0.0)
                {
                    disp[i] = width[i];
                    retval.push_back(disp);
                }
                else if(v1[i] >= width[i] || v2[i] >= width[i] || v3[i] >= width[i])
                {
                    disp[i] = -width[i];
                    retval.push_back(disp);
                }
            }
            return retval;
        };

    std::vector<FaceID> potentially_overlapping_fids;

    for(auto&& fid : this->polygon_->list_face_ids())
    {
        const auto& triangle = this->polygon_->triangle_at(fid);

        Real3 lower(0,0,0), upper(width);
        triangle.bounding_box(width, lower, upper);

        if(shape->test_AABB(lower, upper))
        {
            potentially_overlapping_fids.push_back(fid);
            continue;
        }

        const auto transposes = transpose_direction(triangle);
        for(const auto& transpose : transposes)
        {
            auto vertices = triangle.vertices();
            for(auto& vertex : vertices)
            {
                vertex += transpose;
            }
            const Triangle transposed(vertices);

            transposed.bounding_box(width, lower, upper);
            if(shape->test_AABB(lower, upper))
            {
                potentially_overlapping_fids.push_back(fid);
                break;
            }
        }
    }

    const auto& candidate_faces  = potentially_overlapping_fids;
    const Integer num_candidates = candidate_faces.size();
    if(num_candidates == 0)
    {
        throw std::invalid_argument("The shape does not overlap with polygon.");
    }

    const auto info = this->get_molecule_info(sp);
    for(Integer i=0; i<num; ++i)
    {
        bool particle_inserted = false;
        while(!particle_inserted)
        {
            const auto& face_id = candidate_faces.at(
                    this->rng_->uniform_int(0, num_candidates-1)
                );
            const Real3 pos = this->polygon_->draw_position_on_face(
                    this->rng_, face_id
                );

            if(shape->is_inside(pos))
            {
                const Particle p(sp, pos, info.radius, info.D);
                std::tie(std::ignore, particle_inserted) =
                    this->new_particle(p, face_id);
            }
            // otherwise, particle_inserted is kept false.
        }
    }
    return;
}


std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld::list_particles_within_radius(
        const std::pair<Real3, FaceID>& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>> retval;

    // look particles on the same face
    for(const auto& pid : this->list_particleIDs(pos.second))
    {
        const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
        const Real dist = length(pos.first - pp.second.position()) -
                          pp.second.radius();
        if(dist < radius)
        {
            retval.push_back(std::make_pair(pp, dist));
        }
    }

    // look particles around
    for(const auto& fid : polygon_->neighbor_faces_of(pos.second))
    {
        for(const auto& pid : this->list_particleIDs(fid))
        {
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
            const Real dist = ecell4::polygon::distance(*polygon_,
                pos, std::make_pair(pp.second.position(), get_face_id(pp.first))
                ) - pp.second.radius();
            if(dist < radius)
            {
                retval.push_back(std::make_pair(pp, dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld::list_particles_within_radius(
        const std::pair<Real3, FaceID>& pos, const Real& radius,
        const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real>> retval;

    // look particles on the same face
    for(const auto& pid : this->list_particleIDs(pos.second))
    {
        if(pid == ignore)
        {
            continue;
        }
        const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
        const Real dist = length(pos.first - pp.second.position()) -
                          pp.second.radius();

        if(dist < radius)
        {
            retval.push_back(std::make_pair(pp, dist));
        }
    }

    for(const auto& fid : polygon_->neighbor_faces_of(pos.second))
    {
        for(const auto& pid : this->list_particleIDs(fid))
        {
            if(pid == ignore)
            {
                continue;
            }
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
            const Real dist = ecell4::polygon::distance(*polygon_,
                pos, std::make_pair(pp.second.position(), get_face_id(pp.first))
                ) - pp.second.radius();
            if(dist < radius)
            {
                retval.push_back(std::make_pair(pp, dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}

std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
SGFRDWorld::list_particles_within_radius(
        const std::pair<Real3, FaceID>& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    {// same face
        for(const auto& pid : this->list_particleIDs(pos.second))
        {
            if(pid == ignore1 || pid == ignore2)
            {
                continue;
            }
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius)
            {
                retval.push_back(std::make_pair(pp, dist));
            }
        }
    }

    for(const auto& fid : polygon_->neighbor_faces_of(pos.second))
    {
        for(const auto& pid : this->list_particleIDs(fid))
        {
            if(pid == ignore1 || pid == ignore2)
            {
                continue;
            }
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(pid);
            const Real dist = ecell4::polygon::distance(*polygon_,
                pos, std::make_pair(pp.second.position(), get_face_id(pp.first))
                ) - pp.second.radius();
            if(dist < radius)
            {
                retval.push_back(std::make_pair(pp, dist));
            }
        }
    }
    std::sort(retval.begin(), retval.end(),
              ecell4::utils::pair_second_element_comparator<
                  std::pair<ParticleID, Particle>, Real>());
    return retval;
}


bool SGFRDWorld::check_no_overlap(
        const std::pair<Real3, FaceID>& pos, const Real& radius) const
{
    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(std::vector<ParticleID>::const_iterator
            i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) {return false;}
        }
    }

    std::vector<FaceID> const& neighbors =
        polygon_->neighbor_faces_of(pos.second);
    for(std::vector<FaceID>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = registrator_.elements_over(*iter);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = ecell4::polygon::distance(*polygon_, pos,
                std::make_pair(pp.second.position(), get_face_id(pp.first))) -
                pp.second.radius();
            if(dist < radius) {return false;}
        }
    }
    return true;
}

bool SGFRDWorld::check_no_overlap(
        const std::pair<Real3, FaceID>& pos, const Real& radius,
        const ParticleID& ignore) const
{
    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(std::vector<ParticleID>::const_iterator
            i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore) {continue;}
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) {return false;}
        }
    }

    std::vector<FaceID> const& neighbors =
        polygon_->neighbor_faces_of(pos.second);
    for(std::vector<FaceID>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = this->list_particleIDs(*iter);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore) {continue;}
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = ecell4::polygon::distance(*polygon_, pos,
                std::make_pair(pp.second.position(), get_face_id(pp.first))) -
                pp.second.radius();
            if(dist < radius) {return false;}
        }
    }
    return true;
}

bool SGFRDWorld::check_no_overlap(
        const std::pair<Real3, FaceID>& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const
{
    {// same face
        const std::vector<ParticleID>& ids = this->list_particleIDs(pos.second);
        for(std::vector<ParticleID>::const_iterator
            i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore1 || *i == ignore2) {continue;}
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = length(pos.first - pp.second.position()) -
                              pp.second.radius();
            if(dist < radius) {return false;}
        }
    }

    std::vector<FaceID> const& neighbors =
        polygon_->neighbor_faces_of(pos.second);
    for(std::vector<FaceID>::const_iterator
        iter = neighbors.begin(); iter != neighbors.end(); ++iter)
    {
        const std::vector<ParticleID>& ids = this->list_particleIDs(*iter);
        for(std::vector<ParticleID>::const_iterator
                i(ids.begin()), e(ids.end()); i != e; ++i)
        {
            if(*i == ignore1 || *i == ignore2) {continue;}
            const std::pair<ParticleID, Particle> pp = ps_->get_particle(*i);
            const Real dist = ecell4::polygon::distance(*polygon_, pos,
                std::make_pair(pp.second.position(), get_face_id(pp.first))) -
                pp.second.radius();
            if(dist < radius) {return false;}
        }
    }
    return true;
}

}// sgfrd
}// ecell4
