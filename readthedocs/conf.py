import sys
import os

extensions = ['sphinx.ext.pngmath', 'sphinx.ext.autodoc', 'sphinx.ext.autosummary', 'numpydoc', 'IPython.sphinxext.ipython_console_highlighting', 'IPython.sphinxext.ipython_directive']

templates_path = ['_templates']

master_doc = 'index'

project = u'E-Cell'
copyright = u'2015-, E-Cell project'
author = u'Kazunari Kaizu'

version = '4.1.2'
release = '4.1.2'

language = None

exclude_patterns = ['_build']

pygments_style = 'sphinx'

import sphinx_rtd_theme
html_theme = "sphinx_rtd_theme"
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_static_path = ['_static']

html_show_sourcelink = True

htmlhelp_basename = 'Testdoc'

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',

# Latex figure (float) alignment
#'figure_align': 'htbp',
}

latex_documents = [
  (master_doc, 'Test.tex', u'Test Documentation',
   u'Test', 'manual'),
]

man_pages = [
    (master_doc, 'test', u'Test Documentation',
     [author], 1)
]

texinfo_documents = [
  (master_doc, 'Test', u'Test Documentation',
   author, 'Test', 'One line description of project.',
   'Miscellaneous'),
]

from recommonmark.parser import CommonMarkParser

# The suffix of source filenames.
source_suffix = ['.rst', '.md']

source_parsers = {
	'.md': CommonMarkParser,
}
