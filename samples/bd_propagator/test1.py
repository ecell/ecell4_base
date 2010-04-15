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

w = _gfrd.World(1., 10)
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
    propagator = _gfrd.BDPropagator(w, t, nrw, myrandom.rng, 1e-3, particle_id_list)
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

