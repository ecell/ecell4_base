  <span class="target" id="module-ecell4.util.viz"></span><p>ecell4.util.viz: Visualizer of particles based on D3.js, THREE.js
and Elegans.</p>
<dl class="class">
<dt id="ecell4.util.viz.ColorScale">
<em class="property">class </em><code class="descclassname">ecell4.util.viz.</code><code class="descname">ColorScale</code><span class="sig-paren">(</span><em>config={}</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.viz.ColorScale" title="Permalink to this definition">¶</a></dt>
<dd><p>Color scale for species.</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.util.viz.ColorScale.get_color" title="ecell4.util.viz.ColorScale.get_color"><code class="xref py py-obj docutils literal"><span class="pre">get_color</span></code></a>(name)</td>
<td>Get color unique to the recieved name</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.util.viz.ColorScale.get_config" title="ecell4.util.viz.ColorScale.get_config"><code class="xref py py-obj docutils literal"><span class="pre">get_config</span></code></a>()</td>
<td>Get an instance of dic as the config of colors.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.util.viz.ColorScale.get_color">
<code class="descname">get_color</code><span class="sig-paren">(</span><em>name</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.viz.ColorScale.get_color" title="Permalink to this definition">¶</a></dt>
<dd><p>Get color unique to the recieved name</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>name</strong> : string</p>
<blockquote class="last">
<div><p>This method returns one color unique to this parameter.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.util.viz.ColorScale.get_config">
<code class="descname">get_config</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.viz.ColorScale.get_config" title="Permalink to this definition">¶</a></dt>
<dd><p>Get an instance of dic as the config of colors.</p>
</dd></dl>

</dd></dl>

<dl class="function">
<dt id="ecell4.util.viz.generate_html">
<code class="descclassname">ecell4.util.viz.</code><code class="descname">generate_html</code><span class="sig-paren">(</span><em>keywords</em>, <em>tmpl_path</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.viz.generate_html" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate static html file from JSON model and its own id.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>model</strong> : dict</p>
<blockquote>
<div><p>JSON model from which ecell4.viz generates a plot.</p>
</div></blockquote>
<p><strong>model_id</strong> : string</p>
<blockquote>
<div><p>Unique id for the plot.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">html :</p>
<blockquote class="last">
<div><p>A HTML object</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.util.viz.plot_dense_array">
<code class="descclassname">ecell4.util.viz.</code><code class="descname">plot_dense_array</code><span class="sig-paren">(</span><em>arr, length=256, ranges=None, colors=['#a6cee3', '#fb9a99'], save_image=False, grid=False</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.viz.plot_dense_array" title="Permalink to this definition">¶</a></dt>
<dd><p>Volume renderer</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arr</strong> : list of numpy.array</p>
<blockquote>
<div><p>i.e. [array([[1,2,3], [2,3,4]]), array([[1,2,3]])]</p>
</div></blockquote>
<p><strong>ranges</strong> : list of tuple</p>
<blockquote>
<div><p>ranges for x, y, and z axis
i.e. [(-100, 100), (-100, 100), (-100, 100)]</p>
</div></blockquote>
<p><strong>colors</strong> : list of string</p>
<blockquote>
<div><p>colors for species</p>
</div></blockquote>
<p><strong>length</strong> : int</p>
<blockquote class="last">
<div><p>length of the texture
256 or 64</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.util.viz.plot_movie">
<code class="descclassname">ecell4.util.viz.</code><code class="descname">plot_movie</code><span class="sig-paren">(</span><em>worlds</em>, <em>radius=None</em>, <em>width=500</em>, <em>height=500</em>, <em>config={}</em>, <em>grid=False</em>, <em>species_list=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.viz.plot_movie" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate a movie from received instances of World and show them
on IPython notebook.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>worlds</strong> : list of World</p>
<blockquote>
<div><p>Worlds to render.</p>
</div></blockquote>
<p><strong>radius</strong> : float, default None</p>
<blockquote>
<div><p>If this value is set, all particles in the world will be rendered
as if their radius are the same.</p>
</div></blockquote>
<p><strong>width</strong> : float, default 500</p>
<blockquote>
<div><p>Width of the plotting area.</p>
</div></blockquote>
<p><strong>height</strong> : float, default 500</p>
<blockquote>
<div><p>Height of the plotting area.</p>
</div></blockquote>
<p><strong>config</strong> : dict, default {}</p>
<blockquote>
<div><p>Dict for configure default colors. Its values are colors unique
to each speices.
Colors included in config dict will never be used for other speices.</p>
</div></blockquote>
<p><strong>species_list</strong> : array of string, default None</p>
<blockquote>
<div><p>If set, plot_movie will not search the list of species</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>cfg</strong> : dict</p>
<blockquote class="last">
<div><p>The config data used in this plot.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.util.viz.plot_number_observer">
<code class="descclassname">ecell4.util.viz.</code><code class="descname">plot_number_observer</code><span class="sig-paren">(</span><em>*args</em>, <em>**kwargs</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.viz.plot_number_observer" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate a plot from NumberObservers and show it on IPython notebook
with matplotlib.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>obs</strong> : NumberObserver (e.g. FixedIntervalNumberObserver)</p>
<p><strong>fmt</strong> : str, optional</p>
<p><strong>opt</strong> : dict, optional</p>
<blockquote class="last">
<div><p>matplotlib plot options.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Examples</p>
<div class="highlight-python"><div class="highlight"><pre><span class="gp">&gt;&gt;&gt; </span><span class="n">plot_number_observer</span><span class="p">(</span><span class="n">obs1</span><span class="p">)</span>
<span class="gp">&gt;&gt;&gt; </span><span class="n">plot_number_observer</span><span class="p">(</span><span class="n">obs1</span><span class="p">,</span> <span class="s">&#39;o&#39;</span><span class="p">)</span>
<span class="gp">&gt;&gt;&gt; </span><span class="n">plot_number_observer</span><span class="p">(</span><span class="n">obs1</span><span class="p">,</span> <span class="n">obs2</span><span class="p">,</span> <span class="n">obs3</span><span class="p">,</span> <span class="p">{</span><span class="s">&#39;linewidth&#39;</span><span class="p">:</span> <span class="mi">2</span><span class="p">})</span>
<span class="gp">&gt;&gt;&gt; </span><span class="n">plot_number_observer</span><span class="p">(</span><span class="n">obs1</span><span class="p">,</span> <span class="s">&#39;k-&#39;</span><span class="p">,</span> <span class="n">obs2</span><span class="p">,</span> <span class="s">&#39;k--&#39;</span><span class="p">)</span>
</pre></div>
</div>
</dd></dl>

<dl class="function">
<dt id="ecell4.util.viz.plot_number_observer_with_nya">
<code class="descclassname">ecell4.util.viz.</code><code class="descname">plot_number_observer_with_nya</code><span class="sig-paren">(</span><em>obs</em>, <em>config={}</em>, <em>width=600</em>, <em>height=400</em>, <em>x=None</em>, <em>y=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.viz.plot_number_observer_with_nya" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate a plot from NumberObservers and show it on IPython notebook
with nyaplot.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>obs</strong> : NumberObserver (e.g. FixedIntervalNumberObserver)</p>
<p><strong>config</strong> : dict, optional</p>
<blockquote>
<div><p>A config data for coloring.</p>
</div></blockquote>
<p><strong>width</strong> : int, optional</p>
<p><strong>height</strong> : int, optional</p>
<p><strong>x</strong> : str, optional</p>
<blockquote>
<div><p>A serial for x-axis. If None, x-axis corresponds time.</p>
</div></blockquote>
<p><strong>y</strong> : str or list of str</p>
<blockquote>
<div><p>Serials for y axis.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>cfg</strong> : dict</p>
<blockquote class="last">
<div><p>The config data used in this plot.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.util.viz.plot_trajectory">
<code class="descclassname">ecell4.util.viz.</code><code class="descname">plot_trajectory</code><span class="sig-paren">(</span><em>obs</em>, <em>width=500</em>, <em>height=500</em>, <em>config={}</em>, <em>grid=True</em>, <em>wireframe=False</em>, <em>max_count=10</em>, <em>save_image=False</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.viz.plot_trajectory" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate a plot from received instance of TrajectoryObserver and show it
on IPython notebook.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>obs</strong> : TrajectoryObserver</p>
<blockquote>
<div><p>TrajectoryObserver to render.</p>
</div></blockquote>
<p><strong>width</strong> : float, default 500</p>
<blockquote>
<div><p>Width of the plotting area.</p>
</div></blockquote>
<p><strong>height</strong> : float, default 500</p>
<blockquote>
<div><p>Height of the plotting area.</p>
</div></blockquote>
<p><strong>config</strong> : dict, default {}</p>
<blockquote>
<div><p>Dict for configure default colors. Its values are colors unique
to each particle.
Colors included in config dict will never be used for other particles.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>cfg</strong> : dict</p>
<blockquote class="last">
<div><p>The config data used in this plot.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.util.viz.plot_world">
<code class="descclassname">ecell4.util.viz.</code><code class="descname">plot_world</code><span class="sig-paren">(</span><em>world</em>, <em>radius=None</em>, <em>width=500</em>, <em>height=500</em>, <em>config={}</em>, <em>grid=True</em>, <em>save_image=False</em>, <em>wireframe=False</em>, <em>species_list=None</em>, <em>debug=None</em>, <em>max_count=1000</em>, <em>predicator=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.viz.plot_world" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate a plot from received instance of World and show it on IPython notebook.
This method returns the instance of dict that indicates color setting
for each speices. You can use the dict as the parameter of plot_world,
in order to use the same colors in another plot.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>world</strong> : World</p>
<blockquote>
<div><p>World to render.</p>
</div></blockquote>
<p><strong>radius</strong> : float, default None</p>
<blockquote>
<div><p>If this value is set, all particles in the world will be rendered
as if their radius are the same.</p>
</div></blockquote>
<p><strong>width</strong> : float, default 500</p>
<blockquote>
<div><p>Width of the plotting area.</p>
</div></blockquote>
<p><strong>height</strong> : float, default 500</p>
<blockquote>
<div><p>Height of the plotting area.</p>
</div></blockquote>
<p><strong>config</strong> : dict, default {}</p>
<blockquote>
<div><p>Dict for configure default colors. Its values are colors unique
to each speices.
Colors included in config dict will never be used for other speices.</p>
</div></blockquote>
<p><strong>species_list</strong> : array of string, default None</p>
<blockquote>
<div><p>If set, plot_world will not search the list of species.</p>
</div></blockquote>
<p><strong>debug</strong> : array of dict, default []</p>
<blockquote>
<div><p><strong>* EXPERIMENTAL IMPRIMENTATION *</strong>
Example:
&gt;&gt; [{&#8216;type&#8217;: &#8216;box&#8217;, &#8216;x&#8217;: 10, &#8216;y&#8217;: 10, &#8216;z&#8217;: 10, &#8216;options&#8217;: {&#8216;width&#8217;: 1, &#8216;height&#8217;: 1}}]
type: &#8216;box&#8217;, &#8216;plane&#8217;, &#8216;sphere&#8217;, and &#8216;cylinder&#8217;
x, y, z: float
options:</p>
<blockquote>
<div><p>box: width, height, depth
plane: width, height
sphere: radius
cylinder: radius, height</p>
</div></blockquote>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>cfg</strong> : dict</p>
<blockquote class="last">
<div><p>The config data used in this plot.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>
