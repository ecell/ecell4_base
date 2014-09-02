import uuid
import json
import os

def init_cyjs():
    from IPython.core.display import display, HTML
    path = os.path.abspath(os.path.dirname(__file__)) + '/templates/init_cyjs.js'
    html = open(path).read()
    return display(HTML(html))

def plot_species(species):
    import json
    from jinja2 import Template
    from IPython.core.display import display, HTML

    nodes = [{ 'data': { 'id': unit.name() } } for unit in species.units()]
    print nodes

    path = os.path.abspath(os.path.dirname(__file__)) + '/templates/template.html'
    template = Template(open(path).read())
    html = template.render(nodes=json.dumps(nodes), edges = None, uuid="cy" + str(uuid.uuid4()))
    return display(HTML(html))
