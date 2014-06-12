class ColorScale:
    def from_config(self, config):
        self.config = config
        self.colors = ["#a6cee3","#1f78b4","#b2df8a","#33a02c","#e31a1c", "#8dd3c7", "#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f"]     
        self.buffer = self.colors[:]
        for color in self.config.values():
            if color in self.buffer:
                self.buffer.remove(color)
    
    def get_color(self, type):
        if self.config.get(type) is None:
            if len(self.buffer) == 0:
                self.buffer = self.colors[:]
            color = self.buffer.pop(0)
            self.config[type] = color
            return color
        else :
            return self.config.get(type)

    def to_dict(self):
        return self.config

def init_ipynb():
    import os
    from IPython.core.display import display, HTML

    path = os.path.abspath(os.path.dirname(__file__)) + '/templates/init_ipynb.js'
    html = open(path).read()
    return display(HTML(html))

def plot_world(world, radius=None, width=500, height=500, config={}):
    import uuid
    options = {'width':width, 'height':height}
    color_scale = ColorScale()
    color_scale.from_config(config)

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

        radiuses = [particle[1] for particle in species[k]]
        edge_width = world.edge_lengths()[0]

        if radius is None:
            size = 30/edge_width * radiuses[0]
        else:
            size = 30/edge_width * radius

        plot = {}
        plot['data'] = data
        plot['type'] = "Particles"
        plot['options'] = {'color':color_scale.get_color(k), 'name':k, 'size':size}
        i += 1
        plots.append(plot)

    model = {
        'plots':plots,
        'options':options
    };
    model_id = "\"viz" +  str(uuid.uuid4()) + "\"";
    plot_model(model, model_id)
    return color_scale.to_dict()

def plot_model(model, model_id):
    from IPython.core.display import display, HTML
    import json
    import os
    from jinja2 import Template

    path = os.path.abspath(os.path.dirname(__file__)) + '/templates/particles.tmpl'
    template = Template(open(path).read())
    html = template.render(model=json.dumps(model), model_id = model_id)
    return display(HTML(html));
