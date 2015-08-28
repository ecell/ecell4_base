"""ecell4.util.viz: Visualizer of particles based on D3.js, THREE.js
and Elegans.
"""

import os
import uuid
import json
import base64
import copy
import random
import types


def __parse_world(
        world, radius=None, species_list=None, max_count=None,
        predicator=None):
    """Private function to parse world. Return infomation about particles
    (name, coordinates and particle size) for each species.
    """
    from ecell4 import Species

    if species_list is None:
        species_list = [
            p.species().serial() for pid, p in world.list_particles()]
        species_list = sorted(
            set(species_list), key=species_list.index)  # XXX: pick unique ones

    species = []
    for name in species_list:
        particles = [
            {'pos': p.position(), 'r': p.radius()}
            for pid, p in world.list_particles(Species(name))
            if predicator is None or predicator(pid, p)]
        # particles = [
        #     {'pos': p.position(), 'r': p.radius()}
        #     for pid, p in world.list_particles()
        #     if (p.species().serial() == name and
        #         (predicator is None or predicator(pid, p)))]

        if len(particles) == 0:
            continue

        if max_count is not None and len(particles) > max_count:
            particles = random.sample(particles, max_count)

        data = {
            'x': [p['pos'][0] for p in particles],
            'y': [p['pos'][1] for p in particles],
            'z': [p['pos'][2] for p in particles]
        }

        # assume that all particles belong to one species have the same radius
        r = max([p['r'] for p in particles]) if radius is None else radius
        size = 30.0 / max(world.edge_lengths()) * r

        species.append({
            'name': name,
            'data': data,
            'size': size
        })

    return species


def __get_range_of_world(world):
    edge_lengths = world.edge_lengths()
    max_length = max(tuple(edge_lengths))

    rangex = [(edge_lengths[0] - max_length) * 0.5,
              (edge_lengths[0] + max_length) * 0.5]
    rangey = [(edge_lengths[1] - max_length) * 0.5,
              (edge_lengths[1] + max_length) * 0.5]
    rangez = [(edge_lengths[2] - max_length) * 0.5,
              (edge_lengths[2] + max_length) * 0.5]

    return {'x': rangex, 'y': rangey, 'z': rangez}


def plot_movie(
        worlds, radius=None, width=500, height=500, config={}, grid=False,
        species_list=None):
    """Generate a movie from received instances of World and show them
    on IPython notebook.

    Parameters
    ----------
    worlds : list of World
        Worlds to render.
    radius : float, default None
        If this value is set, all particles in the world will be rendered
        as if their radius are the same.
    width: float, default 500
        Width of the plotting area.
    height: float, default 500
        Height of the plotting area.
    config: dict, default {}
        Dict for configure default colors. Its values are colors unique
        to each speices.
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
        'space_mode': 'wireframe',
        'grid': grid,
        'range': __get_range_of_world(worlds[0])
    }

    model_id = '"movie' + str(uuid.uuid4()) + '"'
    color_scale = ColorScale(config=config)

    display(HTML(generate_html({
        'model_id': model_id,
        'names': json.dumps(list(data.keys())),
        'data': json.dumps(list(data.values())),
        'colors': json.dumps([color_scale.get_color(name)
                              for name in data.keys()]),
        'sizes': json.dumps([sizes[name] for name in data.keys()]),
        'options': json.dumps(options)
    }, '/templates/movie.tmpl')))

    return color_scale.get_config()


def plot_world(
        world, radius=None, width=500, height=500, config={}, grid=True,
        save_image=False, wireframe=False, species_list=None, debug=None, max_count=1000,
        predicator=None):
    """Generate a plot from received instance of World and show it
    on IPython notebook.
    This method returns the instance of dict that indicates color setting
    for each speices. You can use the dict as the parameter of plot_world,
    in order to use the same colors in another plot.

    Parameters
    ----------
    world : World
        World to render.
    radius : float, default None
        If this value is set, all particles in the world will be rendered
        as if their radius are the same.
    width: float, default 500
        Width of the plotting area.
    height: float, default 500
        Height of the plotting area.
    config: dict, default {}
        Dict for configure default colors. Its values are colors unique
        to each speices.
        Colors included in config dict will never be used for other speices.
    species_list: array of string, default None
        If set, plot_world will not search the list of species.
    debug: array of dict, default []
        *** EXPERIMENTAL IMPRIMENTATION ***
        example:
          [{'type': 'box', 'x': 10, 'y': 10, 'z': 10,
            'options': {'width': 1, 'height': 1}}]
        type: 'box', 'plane', 'sphere', and 'cylinder'
        x, y, z: float
        options:
            box: width, height, depth
            plane: width, height
            sphere: radius
            cylinder: radius, height
    """
    from IPython.core.display import display, HTML

    species = __parse_world(world, radius, species_list, max_count, predicator)
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

    if debug is not None:
        data = {'type': [], 'x': [], 'y': [], 'z': [], 'options': []}
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
            'grid': grid,
            'save_image': save_image
        }
    }

    if wireframe:
        model['options']['space_mode'] = 'wireframe'

    model_id = '"viz' + str(uuid.uuid4()) + '"'
    display(HTML(generate_html(
        {'model': json.dumps(model), 'model_id': model_id},
        '/templates/particles.tmpl')))
    return color_scale.get_config()


def plot_dense_array(arr, length=256, ranges=None, colors=["#a6cee3", "#fb9a99"], save_image=False, grid=False):
    """ Volume renderer
    Parameters
    ----------
    arr : list of numpy.array
        i.e. [array([[1,2,3], [2,3,4]]), array([[1,2,3]])]
    ranges : list of tuple
        ranges for x, y, and z axis
        i.e. [(-100, 100), (-100, 100), (-100, 100)]
    colors: list of string
        colors for species
    length: int
        length of the texture
        256 or 64
    """
    import numpy
    from PIL import Image
    from base64 import b64encode
    from tempfile import TemporaryFile
    from math import sqrt
    from IPython.core.display import display, HTML

    # unfold 3d box into 2d grid
    def unfold(arr, dtype=None):
        dtype = arr.dtype if dtype is None else dtype
        i = sqrt(arr.shape[2])
        f_per_row, f_per_column = i, i
        # single channel (luminance)
        try:
            depth, height, width = arr.shape[:]
            arr = arr.reshape((depth*height, width))
            new_arr = numpy.empty((height*f_per_column, width*f_per_row), dtype=dtype)
        # multi channel (RGB)
        except ValueError:
            depth, height, width, channel = arr.shape
            arr = arr.reshape((depth*height, width, channel))
            new_arr = numpy.empty((height*f_per_column, width*f_per_row, channel), dtype=dtype)
        for h in range(0, int(f_per_column)):
            for w in range(0, int(f_per_row)):
                val = arr[(f_per_row*h+w)*height : (f_per_row*h+w+1)*height]
                new_arr[h*height : (h+1)*height, w*width : (w+1)*width] = val
        return new_arr

    def hist(arr, ranges, length, color):
        # create sample
        hist, bins = numpy.histogramdd(arr, bins=tuple([length]*3), range=tuple(ranges))
        # standardize value
        colors = [int(color[1:][i*2:(i+1)*2], 16) for i in range(0, 3)]
        len1d = reduce(lambda val, memo: memo*val, hist.shape, 1)
        arr = [((val/numpy.max(hist))*(hist.copy())).reshape(len1d) for val in colors]
        # add blue and green
        return numpy.array(arr, dtype=numpy.int8).transpose().reshape(tuple(list(hist.shape) + [3]))
    ranges = ranges if ranges is not None else [(numpy.min(a), numpy.max(a)) for a in numpy.array(arr).reshape((sum(map(lambda a: len(a), arr)), 3)).transpose()]

    hist_arr = [hist(a, ranges, length, colors[i]) for i, a in enumerate(arr)]
    compressed = reduce(lambda p, n: p+n, hist_arr)

    img = Image.fromarray(unfold(compressed), "RGB")
    fp = TemporaryFile("r+b")
    img.save(fp, "PNG")
    fp.seek(0)
    encoded_url = "data:image/png;base64," + b64encode(fp.read())
    
    model = {
        'plots': [{
            'type': 'Volume',
            'data': encoded_url,
            'options': {
                'name': "",
                'width': length,
                'height': length,
                'depth': length,
                'f_per_row': sqrt(length),
                'f_per_column': sqrt(length)
            }
        }],
        'options': {
            'grid': grid,
            'save_image': save_image
        }
    }

    model_id = '"viz' + str(uuid.uuid4()) + '"'
    display(HTML(generate_html(
        {'model': json.dumps(model), 'model_id': model_id},
        '/templates/particles.tmpl')))

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


def plot_trajectory(
        obs, width=500, height=500, config={}, grid=True, wireframe=False,
        max_count=10, save_image=False):
    """Generate a plot from received instance of TrajectoryObserver and show it
    on IPython notebook.

    Parameters
    ----------
    obs : TrajectoryObserver
        TrajectoryObserver to render.
    width: float, default 500
        Width of the plotting area.
    height: float, default 500
        Height of the plotting area.
    config: dict, default {}
        Dict for configure default colors. Its values are colors unique
        to each particle.
        Colors included in config dict will never be used for other particles.
    """
    from IPython.core.display import display, HTML

    color_scale = ColorScale(config=config)
    plots = []

    xmin, xmax, ymin, ymax, zmin, zmax = None, None, None, None, None, None

    data = obs.data()
    if max_count is not None and len(data) > max_count:
        data = random.sample(data, max_count)

    for i, y in enumerate(data):
        xarr, yarr, zarr = [], [], []
        for pos in y:
            xarr.append(pos[0])
            yarr.append(pos[1])
            zarr.append(pos[2])

        if xmin is None:
            if len(y) > 0:
                xmin, xmax = min(xarr), max(xarr)
                ymin, ymax = min(yarr), max(yarr)
                zmin, zmax = min(zarr), max(zarr)
        else:
            xmin, xmax = min([xmin] + xarr), max([xmax] + xarr)
            ymin, ymax = min([ymin] + yarr), max([ymax] + yarr)
            zmin, zmax = min([zmin] + zarr), max([zmax] + zarr)

        name = str(i + 1)
        c = color_scale.get_color(name)
        plots.append({
            'type': 'Line',
            'data': {'x': xarr, 'y': yarr, 'z': zarr},
            'options': {
                'name': name,
                'thickness': 2,  # XXX: 'thikness' doesn't work on Windows
                'colors': [c, c]}
        })

    if xmin is None:
        xmin, xmax, ymin, ymax, zmin, zmax = 0, 1, 0, 1, 0, 1

    max_length = max(xmax - xmin, ymax - ymin, zmax - zmin)
    rangex = [(xmin + xmax - max_length) * 0.5,
              (xmin + xmax + max_length) * 0.5]
    rangey = [(ymin + ymax - max_length) * 0.5,
              (ymin + ymax + max_length) * 0.5]
    rangez = [(zmin + zmax - max_length) * 0.5,
              (zmin + zmax + max_length) * 0.5]

    model = {
        'plots': plots,
        'options': {
            'width': width,
            'height': height,
            'range': {'x': rangex, 'y': rangey, 'z': rangez},
            'autorange': False,
            'grid': grid,
            'save_image': save_image
        }
    }

    if wireframe:
        model['options']['space_mode'] = 'wireframe'

    model_id = '"viz' + str(uuid.uuid4()) + '"'
    display(HTML(generate_html(
        {'model': json.dumps(model), 'model_id': model_id},
        '/templates/particles.tmpl')))
    return color_scale.get_config()


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

    filenames = [
       os.path.join(os.path.abspath(os.path.dirname(__file__)),
                    '/templates/ecelllogo/logo%02d.png' % (i + 1))
       for i in range(15)]
    base64s = [
        base64.b64encode(open(filename, 'rt').read())
        for filename in filenames]
    img_html = ('<img name="ecelllogo" style="position:relative;'
                + ' left:0px;" alt="ecelllogo"'
                + ' src="data:image/png;base64,%s"' % (base64s[0])
                + ' onClick="action();" />')
    h = HTML(template % tuple(base64s + [("<p>%s</p>" % (img_html * x)) * y]))
    display(h)

def plot_number_observer(*args, **kwargs):
    """Generate a plot from NumberObservers and show it on IPython notebook
    with matplotlib.

    Parameters
    ----------
    obs : NumberObserver (e.g. FixedIntervalNumberObserver)
    fmt : str
    opt : dict
        matplotlib plot options.

    Examples
    --------
    >>> plot_number_observer(obs1)
    >>> plot_number_observer(obs1, 'o')
    >>> plot_number_observer(obs1, obs2, obs3, {'linewidth': 2})
    >>> plot_number_observer(obs1, 'k-', obs2, 'k--')
    """
    import matplotlib.pylab as plt
    import numpy
    import collections

    special_keys = ("xlim", "ylim", "xlabel", "ylabel", "with_legend", "x", "y")
    plot_opts = {key: value for key, value in kwargs.items()
                 if key not in special_keys}
    color_cycle = plt.rcParams['axes.color_cycle']

    if "y" in kwargs.keys() and isinstance(kwargs["y"], str):
        kwargs["y"] = (kwargs["y"], )

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if len(args) > 1 and isinstance(args[1], str):
        if len(args) % 2 == 0:
            observers = [(args[i], args[i + 1]) for i in range(0, len(args), 2)]
        else:
            observers = [(args[i], args[i + 1]) for i in range(0, len(args) - 1, 2)]
            observers.append(args[-1], None)
    else:
        observers = [(obs, None) for obs in args]

    color_map = {}
    data, xidx = None, 0
    for obs, fmt in observers:
        if isinstance(obs, types.FunctionType):
            if data is None:
                raise ValueError("A function must be given after an observer.")
            y = [obs(xi) for xi in data[xidx]]
            opts = plot_opts.copy()
            opts["label"] = obs.__name__
            if fmt is None:
                ax.plot(data[xidx], y, **opts)
            else:
                ax.plot(data[xidx], y, fmt, **opts)
            continue

        data = numpy.array(obs.data()).T

        if "x" in kwargs.keys():
            targets = [sp.serial() for sp in obs.targets()]
            if kwargs["x"] not in targets:
                raise ValueError("[{0}] given as 'x' was not found.".fomrat(kwargs["x"]))
            xidx = targets.index(kwargs["x"]) + 1
        else:
            xidx = 0

        if "y" in kwargs.keys():
            targets = [sp.serial() for sp in obs.targets()]
            targets = [(targets.index(serial), serial)
                       for serial in kwargs["y"] if serial in targets]
        else:
            targets = [sp.serial() for sp in obs.targets()]
            targets = list(enumerate(targets))
            targets.sort(key=lambda x: x[1])

        for idx, serial in targets:
            opts = plot_opts.copy()

            label = serial
            if len(label) > 0 and label[0] == '_':
                label = '$\_$' + label[1:]  # XXX: lazy escaping for a special character
            if label not in color_map.keys():
                color_map[label] = color_cycle[len(color_map) % len(color_cycle)]
                opts["label"] = label
            opts["color"] = color_map[label]

            if fmt is None:
                ax.plot(data[xidx], data[idx + 1], **opts)
            else:
                ax.plot(data[xidx], data[idx + 1], fmt, **opts)

    if "with_legend" not in kwargs.keys() or kwargs["with_legend"]:
        ax.legend(*ax.get_legend_handles_labels(), loc="best", shadow=True)
    if "xlabel" in kwargs.keys():
        ax.set_xlabel(kwargs["xlabel"])
    elif "x" in kwargs.keys():
        ax.set_xlabel("The Number of Molecules [{0}]".format(kwargs["x"]))
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
    plt.show()

def plot_number_observer_with_nya(obs, config={}, width=600, height=400, x=None, y=None):
    from IPython.core.display import display, HTML
    import numpy

    config = {}
    color_scale = ColorScale(config=config)

    data1, data2 = [], []
    data = numpy.array(obs.data())

    if x is None:
        xidx = 0
    else:
        tmp = [sp.serial() for sp in obs.targets()]
        if x not in tmp:
            raise ValueError("[{0}] given as 'x' was not found.".fomrat(x))
        xidx = tmp.index(x) + 1

    if y is None:
        targets = [sp.serial() for sp in obs.targets()]
        targets = list(enumerate(targets))
        targets.sort(key=lambda x: x[1])
    else:
        if isinstance(y, str):
            y = (y, )
        targets = [sp.serial() for sp in obs.targets()]
        targets = [(targets.index(serial), serial)
                   for serial in y if serial in targets]

    for line in data:
        tmp = {"x": line[xidx]}
        for i, (idx, serial) in enumerate(targets):
            tmp["y{0}".format(i + 1)] = line[idx + 1]
        data1.append(tmp)
    for i, (idx, serial) in enumerate(targets):
        label = serial
        tmp = {"type": "line", "data": "data1",
               "options": {"x": "x", "y": "y{0}".format(i + 1),
                           "stroke_width": 2, "title": label,
                           "color": color_scale.get_color(label)}}
        data2.append(tmp)

    xmin, xmax = data.T[xidx].min(), data.T[xidx].max()
    yview = data.T.take([idx + 1 for idx, serial in targets], axis=0)
    ymin, ymax = yview.min(), yview.max()

    model = {
        "data": {"data1": data1},
        "panes": [{"type": 'rectangular',
                   "diagrams": data2,
                   "options": {"width": width, "height": height, "xrange": [xmin, xmax],
                               "yrange": [ymin, ymax], "legend": True, "zoom": True}}]}
    model_id = 'viz{0:s}'.format(uuid.uuid4())
    display(HTML(generate_html(
        {'model': json.dumps(model), 'model_id': model_id},
        '/templates/nya.tmpl')))
    return color_scale.get_config()

class ColorScale:
    """Color scale for species.
    """

    COLORS = [
        "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#e31a1c", "#8dd3c7",
        "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69",
        "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"]

    def __init__(self, config={}):
        """Initialize a color scale

        Parameters
        ----------
        config : dict, default {}
            Dict for configure default colors. Its values are colors unique
            to each key. Colors included in config will never be used.
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
