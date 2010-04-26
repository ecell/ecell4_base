import _gfrd
import myrandom
import vtk

m = _gfrd.Model()
S0 = m.new_species_type()
S0['D'] = '.01'
S0['radius'] = '.01'
S0['surface'] = 'default'
S1 = m.new_species_type()
S1['D'] = '.01'
S1['radius'] = '.01'
S1['surface'] = 'default'
S2 = m.new_species_type()
S2['D'] = '.01'
S2['radius'] = '.01'
S2['surface'] = 'default'

colors = {
    S0.id: (1., 0., 0.),
    S1.id: (0., 1., 0.),
    S2.id: (1., 1., 0.),
    }

rr = _gfrd.ReactionRule((S0, S1), (S2, ))
rr['k'] = '.01'
m.network_rules.add_reaction_rule(rr)

nrw = _gfrd.NetworkRulesWrapper(m.network_rules)

class MyParticleContainer(_gfrd.ParticleContainer):
    def __init__(self, world_size):
        _gfrd.ParticleContainer.__init__(self)
        self.particles = {}
        self.surfaces = {}
        self.species = {}
        self.pidgen = _gfrd.ParticleIDGenerator(0)
        self.world_size = world_size

    def add_surface(self, surface):
        self.surfaces[surface.id] = surface

    def add_species(self, species):
        self.species[species.id] = species

    def get_surface(self, id):
        return self.surfaces[id]

    def get_species(self, id):
        if isinstance(id, _gfrd.SpeciesType):
            id = id.id
        return self.species[id]

    def new_particle(self, species_id, position):
        new_pid = self.pidgen()
        species = self.get_species(species_id)
        retval = (new_pid, _gfrd.Particle(position, species.radius, species.D, species.id))
        self.update_particle(retval)
        return retval

    def update_particle(self, pid_particle_pair):
        self.particles[pid_particle_pair[0]] = pid_particle_pair[1]
        return False

    def remove_particle(self, pid):
        del self.particles[pid]

    def get_particle(self, pid):
        p = self.particles.get(pid, None)
        if p is None:
            raise NotFound
        return pid, p

    def check_overlap(self, sphere, ignores):
        retval = []
        for pp in self.particles.iteritems():
            if pp[0] in ignores:
                continue
            dist = _gfrd.distance(pp[1].position, sphere.position) - pp[1].radius
            if dist < sphere.radius:
                retval.append((pp, dist))
        retval.sort(lambda a, b: cmp(a[1], b[1]))
        return retval

    def distance(self, x, y):
        return _gfrd.distance_cyclic(x, y, self.world_size)

    def apply_boundary(self, x):
        return _gfrd.apply_boundary(x, self.world_size)

    def cyclic_transpose(self, x, y):
        return _gfrd.cyclic_transpose(x, y, self.world_size)

    def __iter__(self):
        return self.particles.iteritems()

    def create_transaction(self):
        return _gfrd.TransactionImpl(self)

w = MyParticleContainer(1.0)
region = _gfrd._CuboidalRegion("default",
        _gfrd.Box((.5, .5, .5), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.), 1., 1., .1))
w.add_surface(region)

for s in [S0, S1, S2]:
    w.add_species(_gfrd.SpeciesInfo(s.id, float(s['D']), float(s['radius']), s['surface']))

for i in xrange(0, 300):
    w.new_particle([S0, S1][i % 2],
                   [myrandom.uniform(), myrandom.uniform(), myrandom.uniform()])

wn = vtk.vtkRenderWindow()
int = wn.MakeRenderWindowInteractor()
int.Initialize()
int.SetRenderWindow(wn)
r = vtk.vtkRenderer()
wn.AddRenderer(r)

actors = {}

def create_actor(pp):
    s = vtk.vtkSphereSource()
    s.SetRadius(pp[1].radius)
    s.SetCenter(pp[1].position)
    m = vtk.vtkPolyDataMapper()
    m.SetInput(s.GetOutput())
    a = vtk.vtkActor()
    a.GetProperty().SetColor(colors[pp[1].sid])
    a.SetMapper(m)
    r.AddActor(a)
    actors[pp[0]] = (s, a)

def update_actor(pp):
    actors[pp[0]][0].SetCenter(pp[1].position)

def remove_actor(pp):
    r.RemoveActor(actors[pp[0]][1])
    del actors[pp[0]]

for pp in w:
    create_actor(pp) 

anim = []

def callback(*arg):
    t = w.create_transaction()
    particle_id_list = [pair[0] for pair in w]
    propagator = _gfrd._BDPropagator(w, t, nrw, myrandom.rng, 1e-3, 100, particle_id_list)
    propagator.propagate_all()
    for pp in t.added_particles:
        create_actor(pp)
        s = vtk.vtkSphereSource()
        s.SetCenter(pp[1].position)
        s.SetRadius(.01)
        m = vtk.vtkPolyDataMapper()
        m.SetInput(s.GetOutput())
        a = vtk.vtkActor()
        a.GetProperty().SetColor((1., 1., 1.))
        a.GetProperty().SetOpacity(.2)
        a.SetMapper(m)
        r.AddActor(a)
        anim.append((1, s, a))

    for pp in t.removed_particles:
        remove_actor(pp)

    for pp in t.modified_particles:
        update_actor(pp)

    l = len(anim)
    j = 0
    while j < l:
        i, s, a = anim[j]
        if i >= 4:
            r.RemoveActor(a)
            del anim[j]
            l -= 1
            continue
        s.SetRadius(0.04 * i)
        a.GetProperty().SetOpacity(.3 - .05 * i)
        anim[j] = (i + 1, s, a)
        j += 1

    wn.Render()

    del t

int.CreateRepeatingTimer(100)
int.AddObserver('TimerEvent', callback, .0)
int.Start()

