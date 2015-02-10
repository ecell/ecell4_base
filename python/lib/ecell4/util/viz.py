"""

ecell4.util.viz: Visualizer of particles based on D3.js, THREE.js and Elegans.

"""
import os
import uuid
import json
import base64
import copy

import random

def init_ipynb():
    """Load all depending JavaScript libraries to IPython notebook.

    """
    from IPython.core.display import display, HTML

    path = os.path.abspath(os.path.dirname(__file__)) + '/templates/init_ipynb.js'
    html = open(path).read()
    return display(HTML(html))

def __parse_world(world, radius=None, species_list=None, max_count=None):
    """Private function to parse world. Return infomation about particles (name, coordinates and particle size) for each species.
    """

    if species_list is None:
        species_list = [p.species().serial() for pid, p in world.list_particles()]
        species_list = sorted(set(species_list), key=species_list.index) # pick unique ones

    species = []
    for name in species_list:
        particles = [{'pos': p.position(), 'r': p.radius()}
            for pid, p in world.list_particles() if p.species().serial() == name]

        if max_count is not None and len(particles) > max_count:
            particles = random.sample(particles, max_count)

        data = {
            'x': [p['pos'][0] for p in particles],
            'y': [p['pos'][1] for p in particles],
            'z': [p['pos'][2] for p in particles]
        }

        # assume that all particles belong to one species have the same radius
        r = max([p['r'] for p in particles]) if radius is None else radius
        size = 30/min(world.edge_lengths()) * r

        species.append({
            'name': name,
            'data': data,
            'size': size
        })

    return species

def __get_range_of_world(world):
    edge_lengths = world.edge_lengths()
    max_length = max(tuple(edge_lengths))
    rangex = [(edge_lengths[0] - max_length) * 0.5, (edge_lengths[0] + max_length) * 0.5]
    rangey = [(edge_lengths[1] - max_length) * 0.5, (edge_lengths[1] + max_length) * 0.5]
    rangez = [(edge_lengths[2] - max_length) * 0.5, (edge_lengths[2] + max_length) * 0.5]
    return {'x': rangex, 'y': rangey, 'z': rangez}

def plot_movie(worlds, radius=None, width=500, height=500, config={}, grid=False, species_list=None):
    """Generate a movie from received instances of World and show them on IPython notebook.

    Parameters
    ----------
    worlds : list of World
        Worlds to render.
    radius : float, default None
        If this value is set, all particles in the world will be rendered as if their radius are the same.
    width: float, default 500
        Width of the plotting area.
    height: float, default 500
        Height of the plotting area.
    config: dict, default {}
        Dict for configure default colors. Its values are colors unique to each speices.
        Colors included in config dict will never be used for other speices.
    species_list: array of string, default None
        If set, plot_movie will not search the list of species
    """
    from IPython.core.display import display, HTML
    from jinja2 import Template

    data = {}
    sizes = {}
    for i, world in enumerate(worlds):
        species = __parse_world(world, radius, species_list)
        for species_info in species:
            if data.get(species_info['name']) is None:
                data[species_info['name']] = []
            data[species_info['name']].append({
                'df': species_info['data'],
                't': i
            })
            sizes[species_info['name']] = species_info['size']

    options = {
        'player': True,
        'autorange': False,
        'space_mode':'wireframe',
        'grid': grid,
        'range': __get_range_of_world(worlds[0])
    }

    model_id = "\"movie" +  str(uuid.uuid4()) + "\"";
    color_scale = ColorScale(config=config)

    display(HTML(generate_html({
        'model_id': model_id,
        'names': json.dumps(data.keys()),
        'data': json.dumps(data.values()),
        'colors': json.dumps([color_scale.get_color(name) for name in data.keys()]),
        'sizes': json.dumps([sizes[name] for name in data.keys()]),
        'options': json.dumps(options)
    },'/templates/movie.tmpl')))

    return color_scale.get_config()

def plot_world(world, radius=None, width=500, height=500, config={}, grid=False, species_list=None, debug=None, max_count=1000):
    """Generate a plot from received instance of World and show it on IPython notebook.
    This method returns the instance of dict that indicates color setting for each speices.
    You can use the dict as the parameter of plot_world, in order to use the same colors in another plot.

    Parameters
    ----------
    world : World
        World to render.
    radius : float, default None
        If this value is set, all particles in the world will be rendered as if their radius are the same.
    width: float, default 500
        Width of the plotting area.
    height: float, default 500
        Height of the plotting area.
    config: dict, default {}
        Dict for configure default colors. Its values are colors unique to each speices.
        Colors included in config dict will never be used for other speices.
    species_list: array of string, default None
        If set, plot_world will not search the list of species.
    debug: array of dict, default []
        *** EXPERIMENTAL IMPRIMENTATION ***
        example:
          [{'type': 'box', 'x': 10, 'y': 10, 'z': 10, options:{ 'width': 1, 'height': 1}}]
        type: 'box', 'plane', 'sphere', and 'cylinder'
        x, y, z: float
        options:
            box: width, height, depth
            plane: width, height
            sphere: radius
            cylinder: radius, height
    """
    from IPython.core.display import display, HTML

    species = __parse_world(world, radius, species_list, max_count)
    color_scale = ColorScale(config=config)
    plots = []

    for species_info in species:
        plots.append({
            'type': 'Particles',
            'data': species_info['data'],
            'options': {
                'name': species_info['name'],
                'color': color_scale.get_color(species_info['name']),
                'size': species_info['size']
            }
        })

    if debug != None:
        data = {'type':[], 'x':[], 'y':[], 'z':[], 'options':[]}
        for obj in debug:
            for k, v in obj.items():
                data[k].append(v)

        plots.append({
            'type': 'DebugObject',
            'data': data,
            'options': {}
        })

    model = {
        'plots': plots,
        'options': {
            'width': width,
            'height': height,
            'range': __get_range_of_world(world),
            'autorange': False,
            'space_mode':'wireframe',
            'grid': grid
        }
    }

    model_id = "\"viz" +  str(uuid.uuid4()) + "\"";
    display(HTML(generate_html({'model': json.dumps(model), 'model_id': model_id}, '/templates/particles.tmpl')))
    return color_scale.get_config()

def generate_html(keywords, tmpl_path):
    """Generate static html file from JSON model and its own id.

    Parameters
    ----------
    model : dict
        JSON model from which ecell4.viz generates a plot.
    model_id : string
        Unique id for the plot.
    """
    from jinja2 import Template

    path = os.path.abspath(os.path.dirname(__file__)) + tmpl_path
    template = Template(open(path).read())
    html = template.render(**keywords)
    return html

def logo(x=1, y=None):
    if not isinstance(x, int):
        x = 1
    else:
        x = min(10, max(1, x))
    if y is None or not isinstance(y, int):
        y = 1
    else:
        y = min(10, max(1, y))

    from IPython.core.display import display, HTML, Javascript

    template = """<script type="text/javascript">
    var y = 0;
    var running = false, stop = true;
    var base64a = ["%s", "%s", "%s", "%s", "%s",
        "%s", "%s", "%s", "%s", "%s",
        "%s", "%s", "%s", "%s", "%s"];
    var maxcnt = base64a.length;
    var timer_id;

    function move() {
        if (running)
        {
            y = (y + 1) %% maxcnt;
            var logos = document.getElementsByName('ecelllogo');
            for (var i = 0; i < logos.length; i++) {
                logos[i].src = "data:image/png;base64," + base64a[y + 1];
            }
            if (stop && y == maxcnt - 1) {
                // clearInterval(id);
                running = false;
                stop = true;
            }
        }
    }

    function action() {
        if (!stop) {
            stop = true;
        }
        else if (!running) {
            running = true;
            stop = false;
            if (timer_id != undefined) {
                clearInterval(timer_id);
            }
            timer_id = setInterval('move();', 120);
        }
    }
    </script>
    %s
    """

    filenames = [os.path.abspath(os.path.dirname(__file__))
        + '/templates/ecelllogo/logo%02d.png' % (i + 1) for i in range(15)]
    base64s = [base64.b64encode(open(filename, 'rt').read())
        for filename in filenames]
    img_html = '<img name="ecelllogo" style="position:relative; left:0px;" alt="ecelllogo" src="data:image/png;base64,%s" onClick="action();" />' % (base64s[0])
    h = HTML(template % tuple(base64s + [("<p>%s</p>" % (img_html * x)) * y]))
    # j = Javascript("running = true; stop = false; id = setInterval('move();', %g);" % interval)
    # display(h, j)
    display(h)

def plot_number_observer(*args, **kwargs):
    import matplotlib.pylab as plt
    import numpy
    import collections

    special_keys = ("xlim", "ylim", "xlabel", "ylabel")
    plot_opts = {key: value for key, value in kwargs.items() if key not in special_keys}

    fig = plt.figure()
    ax = fig.add_subplot(111)

    is_first = True
    color_cycle = plt.rcParams['axes.color_cycle']
    if len(args) != 1 and isinstance(args[1], str):
        for obs, fmt in zip(args[: : 2], args[1: : 2]):
            data = numpy.array(obs.data()).T
            for i, sp in enumerate(obs.targets()):
                if is_first:
                    ax.plot(data[0], data[i + 1], fmt,
                        color=color_cycle[i % len(color_cycle)], label=sp.serial(), **plot_opts)
                else:
                    ax.plot(data[0], data[i + 1], fmt,
                        color=color_cycle[i % len(color_cycle)], **plot_opts)
            is_first = False
    else:
        for obs in args:
            data = numpy.array(obs.data()).T
            for i, sp in enumerate(obs.targets()):
                if is_first:
                    ax.plot(data[0], data[i + 1],
                        color=color_cycle[i % len(color_cycle)], label=sp.serial(), **plot_opts)
                else:
                    ax.plot(data[0], data[i + 1],
                        color=color_cycle[i % len(color_cycle)], **plot_opts)
            is_first = False

    ax.legend(*ax.get_legend_handles_labels(), loc="best", shadow=True)
    if "xlabel" in kwargs.keys():
        ax.set_xlabel(kwargs["xlabel"])
    else:
        ax.set_xlabel("Time")
    if "ylabel" in kwargs.keys():
        ax.set_ylabel(kwargs["ylabel"])
    else:
        ax.set_ylabel("The Number of Molecules")
    if "xlim" in kwargs.keys():
        ax.set_xlim(kwargs["xlim"])
    if "ylim" in kwargs.keys():
        ax.set_ylim(kwargs["ylim"])
    # fig.show()

class ColorScale:
    """Color scale for species.
    """

    COLORS = ["#a6cee3","#1f78b4","#b2df8a","#33a02c","#e31a1c", "#8dd3c7", "#ffffb3",
        "#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#d9d9d9",
        "#bc80bd","#ccebc5","#ffed6f"]

    def __init__(self, config={}):
        """Initialize a color scale

        Parameters
        ----------
        config : dict, default {}
            Dict for configure default colors. Its values are colors unique to each key.
            Colors included in config will never be used.
        """

        self.config = config
        self.buffer = ColorScale.COLORS[:]

        for color in self.config.values():
            if color in self.buffer:
                self.buffer.remove(color)

    def get_color(self, name):
        """Get color unique to the recieved name

        Parameters
        ----------
        name : string
            This method returns one color unique to this parameter.
        """

        if self.config.get(name) is None:
            self.config[name] = self.buffer.pop(0)
            if len(self.buffer) == 0:
                self.buffer = ColorScale.COLORS[:]

        return self.config[name]

    def get_config(self):
        """Get an instance of dic as the config of colors.
        """
        return self.config
