"""

ecell4.util.viz: Visualizer of particles based on D3.js, THREE.js and Elegans.

"""
import os
import uuid
import json
import base64
import copy

def init_ipynb():
    """Load all depending JavaScript libraries to IPython notebook.

    """
    from IPython.core.display import display, HTML

    path = os.path.abspath(os.path.dirname(__file__)) + '/templates/init_ipynb.js'
    html = open(path).read()
    return display(HTML(html))

def parse_world(world, radius=None, config, species_list=None):
    if species_list is None:
        species = [p.species().serial() for pid, p in world.list_particles()]
        species = sorted(set(species), key=species.index) # pick unique ones
    else:
        species = copy.copy(species_list)

    color_scale = ColorScale(config=config)

    info = {'particles': [], 'ranges': {}}

    for name in species:
        particles = [{'pos': p.position(), 'r': p.radius()}
            for pid, p in world.list_particles() if p.species().serial() == name]
        data = {
            'x': [p['pos'][0] for p in particles],
            'y': [p['pos'][1] for p in particles],
            'z': [p['pos'][2] for p in particles]
        }

        color = color_scale.get_color(name)

        # assume that all particles have the same radius
        r = max([p['r'] for p in particles]) if radius is None else radius
        size = 30/min(world.edge_lengths()) * r

        info['particles'].append({
            'name': name,
            'color': color,
            'data': data,
            'size': size
        })

    edge_lengths = world.edge_lengths()
    max_length = max(tuple(edge_lengths))
    rangex = [(edge_lengths[0] - max_length) * 0.5, (edge_lengths[0] + max_length) * 0.5]
    rangey = [(edge_lengths[1] - max_length) * 0.5, (edge_lengths[1] + max_length) * 0.5]
    rangez = [(edge_lengths[2] - max_length) * 0.5, (edge_lengths[2] + max_length) * 0.5]

    info['ranges'] = {'x': rangex, 'y': rangey, 'z': rangez}
    config = color_scale.get_config()

    return info, config


def plot_movie(worlds, radius=None, width=500, height=500, config={}, grid=False, species_list=None):
    from IPython.core.display import display, HTML
    from jinja2 import Template

    # find information in world[0]
    info, config = parse_world(worlds[0], radius, config, species_list)
    ranges = info['ranges']
    data = { species['name']: {'data': [], 'name': species['name'], 'color': species['color']} for species in info['particles']}

    # find information in each worlds
    i=0
    for world in worlds:
        info = parse_world(world, config, species_list)
        for species in info['particles']:
            data[species['name']]['data'].append({
                'df': species['data'],
                't': i
            })
        i += 1

    options = {
        'player': True,
        'autorange': False,
        'space_mode':'wireframe',
        'grid': grid
        'range': ranges
    }

    model_id = "\"movie" +  str(uuid.uuid4()) + "\"";

    display(HTML(generate_html({
        'model_id': model_id,
        'data': json.dumps([d['data'] for d in data.values()]),
        'options': json.dumps(options),
        'colors': json.dumps([d['color'] for d in data.values()]),
        'names': json.dumps([d['name'] for d in data.values()]),
    },'/templates/movie.tmpl')))

    return config

def plot_world(world, radius=None, width=500, height=500, config={}, grid=False, species_list=None):
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
    """
    from IPython.core.display import display, HTML

    info, config = parse_world(world, radius, config, species_list)

    plot = []
    for species in info['particles']:
        plot.append({
            'type': 'Particles',
            'data': species['data'],
            'options': {'name': species['name'], 'color': species['color'], 'size': species['size']}
        })

    model = {
        'plots': plots,
        'options': {
            'width': width,
            'height': height,
            'range': info['ranges'],
            'autorange': False,
            'space_mode':'wireframe',
            'grid': grid
        }
    };

    model_id = "\"viz" +  str(uuid.uuid4()) + "\"";
    display(HTML(generate_html({'model': json.dumps(model), 'model_id': model_id}, '/templates/particles.tmpl')))
    return config

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
