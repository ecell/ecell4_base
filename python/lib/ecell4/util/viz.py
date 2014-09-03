"""

ecell4.util.viz: Visualizer of particles based on D3.js, THREE.js and Elegans.

"""
import os
import uuid
import json
import base64

def init_ipynb():
    """Load all depending JavaScript libraries to IPython notebook.

    """
    from IPython.core.display import display, HTML

    path = os.path.abspath(os.path.dirname(__file__)) + '/templates/init_ipynb.js'
    html = open(path).read()
    return display(HTML(html))

def plot_world(world, radius=None, width=500, height=500, config={}):
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

    species = [p.species().serial() for pid, p in world.list_particles()]
    species = sorted(set(species), key=species.index) # pick unique ones

    color_scale = ColorScale(config=config)

    plots = []
    for name in species:
        particles = [{'pos':p.position(), 'r':p.radius()}
            for pid, p in world.list_particles() if p.species().serial() is name]
        data = {
            'x': [p['pos'][0] for p in particles],
            'y': [p['pos'][1] for p in particles],
            'z': [p['pos'][2] for p in particles]
        }

        # assume that all particles has the same radius
        r = max([p['r'] for p in particles]) if radius is None else radius
        size = 30/min(world.edge_lengths()) * r

        plots.append({
            'type':"Particles",
            'data':data,
            'options':{'name':name, 'color':color_scale.get_color(name), 'size':size}
        })

    # model = {
    #     'plots':plots,
    #     'options':{'width': width, 'height': height}
    # };

    edge_lengths = world.edge_lengths()
    max_length = max(tuple(edge_lengths))
    rangex = [(edge_lengths[0] - max_length) * 0.5, (edge_lengths[0] + max_length) * 0.5]
    rangey = [(edge_lengths[1] - max_length) * 0.5, (edge_lengths[1] + max_length) * 0.5]
    rangez = [(edge_lengths[2] - max_length) * 0.5, (edge_lengths[2] + max_length) * 0.5]
    model = {
        'plots': plots,
        'options': {'width': width, 'height': height,
            'range': {'x': rangex, 'y': rangey, 'z': rangez}, 'autorange': False}
    };

    model_id = "\"viz" +  str(uuid.uuid4()) + "\"";
    display(HTML(generate_html(model, model_id)))
    return color_scale.get_config()

def generate_html(model, model_id):
    """Generate static html file from JSON model and its own id.

    Parameters
    ----------
    model : dict
        JSON model from which ecell4.viz generates a plot.
    model_id : string
        Unique id for the plot.
    """
    from jinja2 import Template

    path = os.path.abspath(os.path.dirname(__file__)) + '/templates/particles.tmpl'
    template = Template(open(path).read())
    html = template.render(model=json.dumps(model), model_id = model_id)
    return html

def logo():
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
    </script>
    <p><img name="ecelllogo" style="position:relative; left:0px;" alt="ecelllogo" src="data:image/png;base64,%s" onClick="if (!stop) stop = true; else if (!running) { running = true; stop = false; if (timer_id != undefined) { clearInterval(timer_id); } timer_id = setInterval('move();', 100); }" /></p>
    """

    filenames = [os.path.abspath(os.path.dirname(__file__))
        + '/templates/ecelllogo/logo%02d.png' % (i + 1) for i in range(15)]
    base64s = [base64.b64encode(open(filename, 'rt').read())
        for filename in filenames]
    h = HTML(template % tuple(base64s + [base64s[0]]))
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
