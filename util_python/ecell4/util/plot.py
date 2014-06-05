def init_ipynb():
    import os
    from IPython.core.display import display, HTML

    path = os.path.abspath(os.path.dirname(__file__)) + '/templates/init_ipynb.js'
    html = open(path).read()
    return display(HTML(html))

def plot_world(world, options={}):
    colors = ["#a6cee3","#1f78b4","#b2df8a","#33a02c","#e31a1c"]
    options = dict(options.items() + {'width':500, 'height':500}.items())

    particles = [(tuple(p.position()), p.radius(), p.species().serial()) for pid, p in world.list_particles()]
    species = {}
    for particle in particles:
        if species.get(particle[2]) is None:
            species[particle[2]] = []
        species[particle[2]].append(particle)

    plots = []
    i = 0
    for k in species.keys():
        data = {}
        data['x'] = [particle[0][0] for particle in species[k]]
        data['y'] = [particle[0][1] for particle in species[k]]
        data['z'] = [particle[0][2] for particle in species[k]]

        plot = {}
        plot['data'] = data
        plot['type'] = "Particles"
        plot['options'] = {'color':colors[i], 'name':k}
        i += 1
        plots.append(plot)

    model = {
        'plots':plots,
        'options':options
    };
    plot_model(model)

def plot_model(model):
    from IPython.core.display import display, HTML
    import json
    import os
    from jinja2 import Template

    path = os.path.abspath(os.path.dirname(__file__)) + '/templates/particles.tmpl'
    template = Template(open(path).read())
    html = template.render(model=json.dumps(model))

    return display(HTML(html));
