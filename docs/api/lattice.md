  <span class="target" id="module-ecell4.lattice"></span><dl class="class">
<dt id="ecell4.lattice.LatticeFactory">
<em class="property">class </em><code class="descclassname">ecell4.lattice.</code><code class="descname">LatticeFactory</code><a class="headerlink" href="#ecell4.lattice.LatticeFactory" title="Permalink to this definition">¶</a></dt>
<dd><p>A factory class creating a LatticeWorld instance and a LatticeSimulator instance.</p>
<p>LatticeFactory(Real voxel_radius=None, GSLRandomNumberGenerator rng=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeFactory.create_simulator" title="ecell4.lattice.LatticeFactory.create_simulator"><code class="xref py py-obj docutils literal"><span class="pre">create_simulator</span></code></a>((arg1,&nbsp;...)</td>
<td>Return a LatticeSimulator instance.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeFactory.create_world" title="ecell4.lattice.LatticeFactory.create_world"><code class="xref py py-obj docutils literal"><span class="pre">create_world</span></code></a>((arg1=None)&nbsp;-&gt;&nbsp;LatticeWorld)</td>
<td>Return a LatticeWorld instance.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.lattice.LatticeFactory.create_simulator">
<code class="descname">create_simulator</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2</em><span class="sig-paren">)</span> &rarr; LatticeSimulator<a class="headerlink" href="#ecell4.lattice.LatticeFactory.create_simulator" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a LatticeSimulator instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : LatticeWorld</p>
<blockquote>
<div><p>A world</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : Model</p>
<blockquote>
<div><p>A simulation model</p>
</div></blockquote>
<p><strong>arg2</strong> : LatticeWorld</p>
<blockquote>
<div><p>A world</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">LatticeSimulator:</p>
<blockquote class="last">
<div><p>The created simulator</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeFactory.create_world">
<code class="descname">create_world</code><span class="sig-paren">(</span><em>arg1=None</em><span class="sig-paren">)</span> &rarr; LatticeWorld<a class="headerlink" href="#ecell4.lattice.LatticeFactory.create_world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a LatticeWorld instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : Real3</p>
<blockquote>
<div><p>The lengths of edges of a LatticeWorld created</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : str</p>
<blockquote>
<div><p>The path of a HDF5 file for LatticeWorld</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">LatticeWorld:</p>
<blockquote class="last">
<div><p>The created world</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.lattice.LatticeSimulator">
<em class="property">class </em><code class="descclassname">ecell4.lattice.</code><code class="descname">LatticeSimulator</code><a class="headerlink" href="#ecell4.lattice.LatticeSimulator" title="Permalink to this definition">¶</a></dt>
<dd><p>A class running the simulation with the lattice algorithm.</p>
<p>LatticeSimulator(m, w)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.dt" title="ecell4.lattice.LatticeSimulator.dt"><code class="xref py py-obj docutils literal"><span class="pre">dt</span></code></a></td>
<td>Return the step interval.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.get_alpha" title="ecell4.lattice.LatticeSimulator.get_alpha"><code class="xref py py-obj docutils literal"><span class="pre">get_alpha</span></code></a></td>
<td>Return the value of alpha.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.initialize" title="ecell4.lattice.LatticeSimulator.initialize"><code class="xref py py-obj docutils literal"><span class="pre">initialize</span></code></a></td>
<td>Initialize the simulator.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.last_reactions" title="ecell4.lattice.LatticeSimulator.last_reactions"><code class="xref py py-obj docutils literal"><span class="pre">last_reactions</span></code></a>(()&nbsp;-&gt;&nbsp;[(ReactionRule,&nbsp;...)</td>
<td>Return reactions occuring at the last step.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.model" title="ecell4.lattice.LatticeSimulator.model"><code class="xref py py-obj docutils literal"><span class="pre">model</span></code></a></td>
<td>Return the model bound.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.next_time" title="ecell4.lattice.LatticeSimulator.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the scheduled time for the next step.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.num_steps" title="ecell4.lattice.LatticeSimulator.num_steps"><code class="xref py py-obj docutils literal"><span class="pre">num_steps</span></code></a></td>
<td>Return the number of steps.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.run" title="ecell4.lattice.LatticeSimulator.run"><code class="xref py py-obj docutils literal"><span class="pre">run</span></code></a>(duration,&nbsp;observers)</td>
<td>Run the simulation.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.set_alpha" title="ecell4.lattice.LatticeSimulator.set_alpha"><code class="xref py py-obj docutils literal"><span class="pre">set_alpha</span></code></a>(alpha)</td>
<td>Set the value of alpha.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.set_dt" title="ecell4.lattice.LatticeSimulator.set_dt"><code class="xref py py-obj docutils literal"><span class="pre">set_dt</span></code></a>(dt)</td>
<td>Set a step interval.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.set_t" title="ecell4.lattice.LatticeSimulator.set_t"><code class="xref py py-obj docutils literal"><span class="pre">set_t</span></code></a>(t)</td>
<td>Set the current time.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.step" title="ecell4.lattice.LatticeSimulator.step"><code class="xref py py-obj docutils literal"><span class="pre">step</span></code></a>((upto=None)&nbsp;-&gt;&nbsp;bool)</td>
<td>Step the simulation.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.t" title="ecell4.lattice.LatticeSimulator.t"><code class="xref py py-obj docutils literal"><span class="pre">t</span></code></a></td>
<td>Return the time.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeSimulator.world" title="ecell4.lattice.LatticeSimulator.world"><code class="xref py py-obj docutils literal"><span class="pre">world</span></code></a></td>
<td>Return the world bound.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.dt">
<code class="descname">dt</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.dt" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the step interval.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.get_alpha">
<code class="descname">get_alpha</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.get_alpha" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the value of alpha.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.initialize">
<code class="descname">initialize</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.initialize" title="Permalink to this definition">¶</a></dt>
<dd><p>Initialize the simulator.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.last_reactions">
<code class="descname">last_reactions</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [(ReactionRule, ReactionInfo)]<a class="headerlink" href="#ecell4.lattice.LatticeSimulator.last_reactions" title="Permalink to this definition">¶</a></dt>
<dd><p>Return reactions occuring at the last step.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>The list of reaction rules and infos.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.model">
<code class="descname">model</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.model" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the model bound.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the scheduled time for the next step.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.run">
<code class="descname">run</code><span class="sig-paren">(</span><em>duration</em>, <em>observers</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.run" title="Permalink to this definition">¶</a></dt>
<dd><p>Run the simulation.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>duration</strong> : Real</p>
<blockquote>
<div><p>A duration for running a simulation.
A simulation is expected to be stopped at t() + duration.</p>
</div></blockquote>
<p><strong>observers</strong> : list of Obeservers, optional</p>
<blockquote class="last">
<div><p>observers</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.set_alpha">
<code class="descname">set_alpha</code><span class="sig-paren">(</span><em>alpha</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.set_alpha" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the value of alpha.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>alpha</strong> : Real</p>
<blockquote class="last">
<div><p>The value of alpha</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.set_dt">
<code class="descname">set_dt</code><span class="sig-paren">(</span><em>dt</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.set_dt" title="Permalink to this definition">¶</a></dt>
<dd><p>Set a step interval.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>dt</strong> : Real</p>
<blockquote class="last">
<div><p>A step interval</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.set_t" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the current time.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>t</strong> : Real</p>
<blockquote class="last">
<div><p>A current time.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.step">
<code class="descname">step</code><span class="sig-paren">(</span><em>upto=None</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.lattice.LatticeSimulator.step" title="Permalink to this definition">¶</a></dt>
<dd><p>Step the simulation.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>upto</strong> : Real, optional</p>
<blockquote>
<div><p>The time which to step the simulation up to</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>True if the simulation did not reach the given time.
When upto is not given, nothing will be returned.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeSimulator.world">
<code class="descname">world</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeSimulator.world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the world bound.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.lattice.LatticeWorld">
<em class="property">class </em><code class="descclassname">ecell4.lattice.</code><code class="descname">LatticeWorld</code><a class="headerlink" href="#ecell4.lattice.LatticeWorld" title="Permalink to this definition">¶</a></dt>
<dd><p>A class containing the properties of the lattice world.</p>
<p>LatticeWorld(edge_lengths=None, voxel_radius=None, GSLRandomNumberGenerator rng=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.actual_lengths" title="ecell4.lattice.LatticeWorld.actual_lengths"><code class="xref py py-obj docutils literal"><span class="pre">actual_lengths</span></code></a></td>
<td>Return the actual edge lengths of the world.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.add_molecules" title="ecell4.lattice.LatticeWorld.add_molecules"><code class="xref py py-obj docutils literal"><span class="pre">add_molecules</span></code></a>(sp,&nbsp;num[,&nbsp;shape])</td>
<td>Add some molecules.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.add_structure" title="ecell4.lattice.LatticeWorld.add_structure"><code class="xref py py-obj docutils literal"><span class="pre">add_structure</span></code></a>(sp,&nbsp;shape)</td>
<td>Add a structure.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.as_base" title="ecell4.lattice.LatticeWorld.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Return self as a base class.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.bind_to" title="ecell4.lattice.LatticeWorld.bind_to"><code class="xref py py-obj docutils literal"><span class="pre">bind_to</span></code></a>(m)</td>
<td>Bind a model to the world</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.col_size" title="ecell4.lattice.LatticeWorld.col_size"><code class="xref py py-obj docutils literal"><span class="pre">col_size</span></code></a></td>
<td>Return the size of the column of the world.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.coord2global" title="ecell4.lattice.LatticeWorld.coord2global"><code class="xref py py-obj docutils literal"><span class="pre">coord2global</span></code></a>((coord)&nbsp;-&gt;&nbsp;Integer3)</td>
<td>Transform a coordinate to a global coordinate.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.coord2private" title="ecell4.lattice.LatticeWorld.coord2private"><code class="xref py py-obj docutils literal"><span class="pre">coord2private</span></code></a></td>
<td>Transform a coordinate to a private coordinate.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.coordinate2position" title="ecell4.lattice.LatticeWorld.coordinate2position"><code class="xref py py-obj docutils literal"><span class="pre">coordinate2position</span></code></a>((coord)&nbsp;-&gt;&nbsp;Real3)</td>
<td>Transform a coordinate to a position.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.edge_lengths" title="ecell4.lattice.LatticeWorld.edge_lengths"><code class="xref py py-obj docutils literal"><span class="pre">edge_lengths</span></code></a>(()&nbsp;-&gt;&nbsp;Real3)</td>
<td>Return the edge lengths of the world.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.get_neighbor" title="ecell4.lattice.LatticeWorld.get_neighbor"><code class="xref py py-obj docutils literal"><span class="pre">get_neighbor</span></code></a>((coord,&nbsp;nrand)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the neighbor coordinate of a given coordinate.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.get_neighbor_private" title="ecell4.lattice.LatticeWorld.get_neighbor_private"><code class="xref py py-obj docutils literal"><span class="pre">get_neighbor_private</span></code></a>((coord,&nbsp;nrand)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the neighbor coordinate of a given coordinate in private.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.get_particle" title="ecell4.lattice.LatticeWorld.get_particle"><code class="xref py py-obj docutils literal"><span class="pre">get_particle</span></code></a>(pid)&nbsp;-&gt;&nbsp;(ParticleID,&nbsp;Particle)</td>
<td>Return the particle associated a given ParticleID.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.get_volume" title="ecell4.lattice.LatticeWorld.get_volume"><code class="xref py py-obj docutils literal"><span class="pre">get_volume</span></code></a></td>
<td>Return the actual volume of the world.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.get_voxel" title="ecell4.lattice.LatticeWorld.get_voxel"><code class="xref py py-obj docutils literal"><span class="pre">get_voxel</span></code></a>(arg)&nbsp;-&gt;&nbsp;(ParticleID,&nbsp;Voxle)</td>
<td>Return the voxel having a particle associated with a given ParticleID or coordinate.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.global2coord" title="ecell4.lattice.LatticeWorld.global2coord"><code class="xref py py-obj docutils literal"><span class="pre">global2coord</span></code></a>((g)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Transform a global coordinate to a coordinate.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.global2position" title="ecell4.lattice.LatticeWorld.global2position"><code class="xref py py-obj docutils literal"><span class="pre">global2position</span></code></a>((g)&nbsp;-&gt;&nbsp;Real3)</td>
<td>Transform a global coordinate to a position.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.global2private" title="ecell4.lattice.LatticeWorld.global2private"><code class="xref py py-obj docutils literal"><span class="pre">global2private</span></code></a>((g)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Transform a global coordinate to a private coordinate.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.has_particle" title="ecell4.lattice.LatticeWorld.has_particle"><code class="xref py py-obj docutils literal"><span class="pre">has_particle</span></code></a>((pid)&nbsp;-&gt;&nbsp;bool)</td>
<td>Check if a particle associated with a given particle id exists.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.has_voxel" title="ecell4.lattice.LatticeWorld.has_voxel"><code class="xref py py-obj docutils literal"><span class="pre">has_voxel</span></code></a>((pid)&nbsp;-&gt;&nbsp;bool)</td>
<td>Check if a particle exists.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.layer_size" title="ecell4.lattice.LatticeWorld.layer_size"><code class="xref py py-obj docutils literal"><span class="pre">layer_size</span></code></a></td>
<td>Return the size of layer of the world.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.list_particles" title="ecell4.lattice.LatticeWorld.list_particles"><code class="xref py py-obj docutils literal"><span class="pre">list_particles</span></code></a>((sp)&nbsp;-&gt;&nbsp;[(ParticleID,&nbsp;Particle)])</td>
<td>Return the list of particles.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.list_particles_exact" title="ecell4.lattice.LatticeWorld.list_particles_exact"><code class="xref py py-obj docutils literal"><span class="pre">list_particles_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;[(ParticleID,&nbsp;...)</td>
<td>Return the list of particles of a given species.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.list_voxels" title="ecell4.lattice.LatticeWorld.list_voxels"><code class="xref py py-obj docutils literal"><span class="pre">list_voxels</span></code></a>((sp=None)&nbsp;-&gt;&nbsp;[ParitcleID,&nbsp;Voxel])</td>
<td>Returns the list of voxels.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.list_voxels_exact" title="ecell4.lattice.LatticeWorld.list_voxels_exact"><code class="xref py py-obj docutils literal"><span class="pre">list_voxels_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;[ParitcleID,&nbsp;Voxel])</td>
<td>Returns the list of voxels.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.load" title="ecell4.lattice.LatticeWorld.load"><code class="xref py py-obj docutils literal"><span class="pre">load</span></code></a>(filename)</td>
<td>Load the world from a file.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.new_particle" title="ecell4.lattice.LatticeWorld.new_particle"><code class="xref py py-obj docutils literal"><span class="pre">new_particle</span></code></a>(arg1[,&nbsp;arg2])</td>
<td>Create a new particle.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.new_voxel" title="ecell4.lattice.LatticeWorld.new_voxel"><code class="xref py py-obj docutils literal"><span class="pre">new_voxel</span></code></a>(arg1,&nbsp;arg2)&nbsp;-&gt;&nbsp;(ParticleID,&nbsp;Voxel)</td>
<td>Create a particle.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.num_molecules" title="ecell4.lattice.LatticeWorld.num_molecules"><code class="xref py py-obj docutils literal"><span class="pre">num_molecules</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of molecules.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.num_molecules_exact" title="ecell4.lattice.LatticeWorld.num_molecules_exact"><code class="xref py py-obj docutils literal"><span class="pre">num_molecules_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of molecules of a given species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.num_particles" title="ecell4.lattice.LatticeWorld.num_particles"><code class="xref py py-obj docutils literal"><span class="pre">num_particles</span></code></a>((sp=None)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of particles.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.num_particles_exact" title="ecell4.lattice.LatticeWorld.num_particles_exact"><code class="xref py py-obj docutils literal"><span class="pre">num_particles_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of particles of a given species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.num_voxels" title="ecell4.lattice.LatticeWorld.num_voxels"><code class="xref py py-obj docutils literal"><span class="pre">num_voxels</span></code></a>((sp=None)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of voxels.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.num_voxels_exact" title="ecell4.lattice.LatticeWorld.num_voxels_exact"><code class="xref py py-obj docutils literal"><span class="pre">num_voxels_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of voxels of a given species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.on_structure" title="ecell4.lattice.LatticeWorld.on_structure"><code class="xref py py-obj docutils literal"><span class="pre">on_structure</span></code></a>((sp,&nbsp;coord)&nbsp;-&gt;&nbsp;bool)</td>
<td>Check if the given species would be on the proper structure at the coordinate.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.position2coordinate" title="ecell4.lattice.LatticeWorld.position2coordinate"><code class="xref py py-obj docutils literal"><span class="pre">position2coordinate</span></code></a>((pos)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Transform a position to a coordinate.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.position2global" title="ecell4.lattice.LatticeWorld.position2global"><code class="xref py py-obj docutils literal"><span class="pre">position2global</span></code></a>((pos)&nbsp;-&gt;&nbsp;Integer3)</td>
<td>Transform a position to a global coordinate.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.private2coord" title="ecell4.lattice.LatticeWorld.private2coord"><code class="xref py py-obj docutils literal"><span class="pre">private2coord</span></code></a></td>
<td>Transform a private coordinate to a coordinate.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.private2global" title="ecell4.lattice.LatticeWorld.private2global"><code class="xref py py-obj docutils literal"><span class="pre">private2global</span></code></a>((coord)&nbsp;-&gt;&nbsp;Integer3)</td>
<td>Transform a private coordinate to a global coordinate.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.private2position" title="ecell4.lattice.LatticeWorld.private2position"><code class="xref py py-obj docutils literal"><span class="pre">private2position</span></code></a>((coord)&nbsp;-&gt;&nbsp;Real3)</td>
<td>Transform a private coordinate to a position.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.remove_molecules" title="ecell4.lattice.LatticeWorld.remove_molecules"><code class="xref py py-obj docutils literal"><span class="pre">remove_molecules</span></code></a>(sp,&nbsp;num)</td>
<td>Remove the molecules.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.remove_particle" title="ecell4.lattice.LatticeWorld.remove_particle"><code class="xref py py-obj docutils literal"><span class="pre">remove_particle</span></code></a>(pid)</td>
<td>Remove the particle associated with a given ParticleID.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.remove_voxel" title="ecell4.lattice.LatticeWorld.remove_voxel"><code class="xref py py-obj docutils literal"><span class="pre">remove_voxel</span></code></a>(pid)</td>
<td>Remove the particle associated with a given ParticleID.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.rng" title="ecell4.lattice.LatticeWorld.rng"><code class="xref py py-obj docutils literal"><span class="pre">rng</span></code></a></td>
<td>Return a random number generator object.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.row_size" title="ecell4.lattice.LatticeWorld.row_size"><code class="xref py py-obj docutils literal"><span class="pre">row_size</span></code></a></td>
<td>Return the size of row of the world.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.save" title="ecell4.lattice.LatticeWorld.save"><code class="xref py py-obj docutils literal"><span class="pre">save</span></code></a>(filename)</td>
<td>Save the world to a file.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.set_t" title="ecell4.lattice.LatticeWorld.set_t"><code class="xref py py-obj docutils literal"><span class="pre">set_t</span></code></a>(t)</td>
<td>Set the value of the time of the world.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.shape" title="ecell4.lattice.LatticeWorld.shape"><code class="xref py py-obj docutils literal"><span class="pre">shape</span></code></a>(()&nbsp;-&gt;&nbsp;Integer3)</td>
<td>Return the triplet of sizes of column, row and layer.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.size" title="ecell4.lattice.LatticeWorld.size"><code class="xref py py-obj docutils literal"><span class="pre">size</span></code></a></td>
<td>Return the size of voxels.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.t" title="ecell4.lattice.LatticeWorld.t"><code class="xref py py-obj docutils literal"><span class="pre">t</span></code></a></td>
<td>Return the time of the world.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.update_particle" title="ecell4.lattice.LatticeWorld.update_particle"><code class="xref py py-obj docutils literal"><span class="pre">update_particle</span></code></a>(pid,&nbsp;p)</td>
<td>Update a particle.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.update_voxel" title="ecell4.lattice.LatticeWorld.update_voxel"><code class="xref py py-obj docutils literal"><span class="pre">update_voxel</span></code></a>((pid,&nbsp;v)&nbsp;-&gt;&nbsp;bool)</td>
<td>Update a particle.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.volume" title="ecell4.lattice.LatticeWorld.volume"><code class="xref py py-obj docutils literal"><span class="pre">volume</span></code></a></td>
<td>Return the volume of the world.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.voxel_radius" title="ecell4.lattice.LatticeWorld.voxel_radius"><code class="xref py py-obj docutils literal"><span class="pre">voxel_radius</span></code></a></td>
<td>Return the voxel radius.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.lattice.LatticeWorld.voxel_volume" title="ecell4.lattice.LatticeWorld.voxel_volume"><code class="xref py py-obj docutils literal"><span class="pre">voxel_volume</span></code></a></td>
<td>Return the volume of a voxel.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.actual_lengths">
<code class="descname">actual_lengths</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.actual_lengths" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the actual edge lengths of the world.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real3:</p>
<blockquote class="last">
<div><p>The actual edge lengths of the world</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.add_molecules">
<code class="descname">add_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em>, <em>shape=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.add_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Add some molecules.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species of molecules to add</p>
</div></blockquote>
<p><strong>num</strong> : Integer</p>
<blockquote>
<div><p>The number of molecules to add</p>
</div></blockquote>
<p><strong>shape</strong> : Shape, optional</p>
<blockquote class="last">
<div><p>A shape to add molecules on</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.add_structure">
<code class="descname">add_structure</code><span class="sig-paren">(</span><em>sp</em>, <em>shape</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.add_structure" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a structure.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species suggesting the shape.</p>
</div></blockquote>
<p><strong>shape</strong> : Shape</p>
<blockquote class="last">
<div><p>A shape of the structure.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Return self as a base class. Only for developmental use.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.bind_to">
<code class="descname">bind_to</code><span class="sig-paren">(</span><em>m</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.bind_to" title="Permalink to this definition">¶</a></dt>
<dd><p>Bind a model to the world</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>m</strong> : Model</p>
<blockquote class="last">
<div><p>A model to bind</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.col_size">
<code class="descname">col_size</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.col_size" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the size of the column of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.coord2global">
<code class="descname">coord2global</code><span class="sig-paren">(</span><em>coord</em><span class="sig-paren">)</span> &rarr; Integer3<a class="headerlink" href="#ecell4.lattice.LatticeWorld.coord2global" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a coordinate to a global coordinate.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.coord2private">
<code class="descname">coord2private</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.coord2private" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a coordinate to a private coordinate.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.coordinate2position">
<code class="descname">coordinate2position</code><span class="sig-paren">(</span><em>coord</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.lattice.LatticeWorld.coordinate2position" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a coordinate to a position.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.edge_lengths">
<code class="descname">edge_lengths</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.lattice.LatticeWorld.edge_lengths" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the edge lengths of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.get_neighbor">
<code class="descname">get_neighbor</code><span class="sig-paren">(</span><em>coord</em>, <em>nrand</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.get_neighbor" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the neighbor coordinate of a given coordinate.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>coord</strong> : Integer</p>
<blockquote>
<div><p>A coordinate of a voxel</p>
</div></blockquote>
<p><strong>nrand</strong> : Integer</p>
<blockquote>
<div><p>A key in the range from 0 to 11 to assign a neighbor voxel</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The coordinate of the neighbor voxel</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.get_neighbor_private">
<code class="descname">get_neighbor_private</code><span class="sig-paren">(</span><em>coord</em>, <em>nrand</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.get_neighbor_private" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the neighbor coordinate of a given coordinate in private.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>coord</strong> : Integer</p>
<blockquote>
<div><p>A private coordinate of a voxel</p>
</div></blockquote>
<p><strong>nrand</strong> : Integer</p>
<blockquote>
<div><p>A key in the range from 0 to 11 to assign a neighbor voxel</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The private coordinate of the neighbor voxel</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.get_particle">
<code class="descname">get_particle</code><span class="sig-paren">(</span><em>pid) -&gt; (ParticleID</em>, <em>Particle</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.get_particle" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the particle associated a given ParticleID.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pid</strong> : ParticleID</p>
<blockquote>
<div><p>A id of the particle you want</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">tuple:</p>
<blockquote class="last">
<div><p>A pair of ParticleID and Particle</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.get_volume">
<code class="descname">get_volume</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.get_volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the actual volume of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.get_voxel">
<code class="descname">get_voxel</code><span class="sig-paren">(</span><em>arg) -&gt; (ParticleID</em>, <em>Voxle</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.get_voxel" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the voxel having a particle associated with a given ParticleID
or coordinate.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg</strong> : ParticleID or Integer</p>
<blockquote>
<div><p>An id or coordiante of the particle in the voxel you want</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">tuple:</p>
<blockquote class="last">
<div><p>A pair of ParticleID and Voxel</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.global2coord">
<code class="descname">global2coord</code><span class="sig-paren">(</span><em>g</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.global2coord" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a global coordinate to a coordinate.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>g</strong> : Integer3</p>
<blockquote>
<div><p>A global coordinate</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>A coordinate</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.global2position">
<code class="descname">global2position</code><span class="sig-paren">(</span><em>g</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.lattice.LatticeWorld.global2position" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a global coordinate to a position.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>g</strong> : Integer3</p>
<blockquote>
<div><p>A global coordinate</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real3:</p>
<blockquote class="last">
<div><p>A position</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.global2private">
<code class="descname">global2private</code><span class="sig-paren">(</span><em>g</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.global2private" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a global coordinate to a private coordinate.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>g</strong> : Integer3</p>
<blockquote>
<div><p>A global coordinate</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>A private coordinate</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.has_particle">
<code class="descname">has_particle</code><span class="sig-paren">(</span><em>pid</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.lattice.LatticeWorld.has_particle" title="Permalink to this definition">¶</a></dt>
<dd><p>Check if a particle associated with a given particle id exists.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pid</strong> : ParticleID</p>
<blockquote>
<div><p>A particle id to check</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>if a particle exists, this is true. Otherwise false</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.has_voxel">
<code class="descname">has_voxel</code><span class="sig-paren">(</span><em>pid</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.lattice.LatticeWorld.has_voxel" title="Permalink to this definition">¶</a></dt>
<dd><p>Check if a particle exists.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pid</strong> : ParticleID</p>
<blockquote>
<div><p>A particle id of the particle to check</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>whether a particle associated with a given particle id exists</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.layer_size">
<code class="descname">layer_size</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.layer_size" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the size of layer of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.list_particles">
<code class="descname">list_particles</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.lattice.LatticeWorld.list_particles" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the list of particles.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species, optional</p>
<blockquote>
<div><p>The species of particles to list up
If no species is given, return the whole list of particles.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>The list of particles (of the given species)</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.list_particles_exact">
<code class="descname">list_particles_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.lattice.LatticeWorld.list_particles_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the list of particles of a given species.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>The species of particles to list up</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>The list of particles of a given species</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.list_voxels">
<code class="descname">list_voxels</code><span class="sig-paren">(</span><em>sp=None</em><span class="sig-paren">)</span> &rarr; [ParitcleID, Voxel]<a class="headerlink" href="#ecell4.lattice.LatticeWorld.list_voxels" title="Permalink to this definition">¶</a></dt>
<dd><p>Returns the list of voxels.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species, optional</p>
<blockquote>
<div><p>A species of particles to list up.
If no species is given, return a list of all voxels.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>The list of the pair of ParticleID and Voxel</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.list_voxels_exact">
<code class="descname">list_voxels_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [ParitcleID, Voxel]<a class="headerlink" href="#ecell4.lattice.LatticeWorld.list_voxels_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>Returns the list of voxels.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species, optional</p>
<blockquote>
<div><p>A species of particles to list up.
If no species is given, return a list of all voxels.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>The list of the pair of ParticleID and Voxel</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.load">
<code class="descname">load</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.load" title="Permalink to this definition">¶</a></dt>
<dd><p>Load the world from a file.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>filename</strong> : str</p>
<blockquote class="last">
<div><p>A filename to load from</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.new_particle">
<code class="descname">new_particle</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2=None) -&gt; (ParticleID</em>, <em>Particle</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.new_particle" title="Permalink to this definition">¶</a></dt>
<dd><p>Create a new particle.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : Particle</p>
<blockquote>
<div><p>A particle to be placed.</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : Species</p>
<blockquote>
<div><p>A species of a particle</p>
</div></blockquote>
<p><strong>arg2</strong> : Real3</p>
<blockquote>
<div><p>A coordinate to place a particle</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">tuple:</p>
<blockquote class="last">
<div><p>A pair of ParticleID and Particle of a new particle</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.new_voxel">
<code class="descname">new_voxel</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2) -&gt; (ParticleID</em>, <em>Voxel</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.new_voxel" title="Permalink to this definition">¶</a></dt>
<dd><p>Create a particle.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : Voxel</p>
<blockquote>
<div><p>The information to create</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : Species</p>
<blockquote>
<div><p>The Species of particles to create</p>
</div></blockquote>
<p><strong>arg2</strong> : Integer</p>
<blockquote>
<div><p>The number of particles(voxels)</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">tuple:</p>
<blockquote class="last">
<div><p>A pair of ParticleID and Voxel</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.num_molecules">
<code class="descname">num_molecules</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.num_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of molecules.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species whose molecules you count</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of molecules (of a given species)</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.num_molecules_exact">
<code class="descname">num_molecules_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.num_molecules_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of molecules of a given species.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species whose molecules you count</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of molecules of a given species</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.num_particles">
<code class="descname">num_particles</code><span class="sig-paren">(</span><em>sp=None</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.num_particles" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of particles.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species, optional</p>
<blockquote>
<div><p>The species of particles to count
If no species is given, return the total number of particles.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of particles (of the given species)</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.num_particles_exact">
<code class="descname">num_particles_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.num_particles_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of particles of a given species.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>The species of particles to count</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of particles of a given species</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.num_voxels">
<code class="descname">num_voxels</code><span class="sig-paren">(</span><em>sp=None</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.num_voxels" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of voxels.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species, optional</p>
<blockquote>
<div><p>The species of particles to count</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of voxels (of the given species)</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.num_voxels_exact">
<code class="descname">num_voxels_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.num_voxels_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of voxels of a given species.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>The species of particles to count</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of voxels of a given species</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.on_structure">
<code class="descname">on_structure</code><span class="sig-paren">(</span><em>sp</em>, <em>coord</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.lattice.LatticeWorld.on_structure" title="Permalink to this definition">¶</a></dt>
<dd><p>Check if the given species would be on the proper structure at the coordinate.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species scheduled to be placed</p>
</div></blockquote>
<p><strong>coord</strong> : Integer</p>
<blockquote>
<div><p>A coordinate to be occupied</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>if it is on the proper structure, or not</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.position2coordinate">
<code class="descname">position2coordinate</code><span class="sig-paren">(</span><em>pos</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.lattice.LatticeWorld.position2coordinate" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a position to a coordinate.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>A coordinate</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.position2global">
<code class="descname">position2global</code><span class="sig-paren">(</span><em>pos</em><span class="sig-paren">)</span> &rarr; Integer3<a class="headerlink" href="#ecell4.lattice.LatticeWorld.position2global" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a position to a global coordinate.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer3:</p>
<blockquote class="last">
<div><p>A global coordinate</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.private2coord">
<code class="descname">private2coord</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.private2coord" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a private coordinate to a coordinate.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.private2global">
<code class="descname">private2global</code><span class="sig-paren">(</span><em>coord</em><span class="sig-paren">)</span> &rarr; Integer3<a class="headerlink" href="#ecell4.lattice.LatticeWorld.private2global" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a private coordinate to a global coordinate.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.private2position">
<code class="descname">private2position</code><span class="sig-paren">(</span><em>coord</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.lattice.LatticeWorld.private2position" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a private coordinate to a position.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.remove_molecules">
<code class="descname">remove_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.remove_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove the molecules.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species whose molecules to remove</p>
</div></blockquote>
<p><strong>num</strong> : Integer</p>
<blockquote class="last">
<div><p>A number of molecules to be removed</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.remove_particle">
<code class="descname">remove_particle</code><span class="sig-paren">(</span><em>pid</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.remove_particle" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove the particle associated with a given ParticleID.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pid</strong> : ParticleID</p>
<blockquote class="last">
<div><p>A id of particle to remove</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.remove_voxel">
<code class="descname">remove_voxel</code><span class="sig-paren">(</span><em>pid</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.remove_voxel" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove the particle associated with a given ParticleID.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pid</strong> : ParticleID</p>
<blockquote class="last">
<div><p>A id of particle to remove</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.rng">
<code class="descname">rng</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.rng" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a random number generator object.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.row_size">
<code class="descname">row_size</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.row_size" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the size of row of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.save">
<code class="descname">save</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.save" title="Permalink to this definition">¶</a></dt>
<dd><p>Save the world to a file.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>filename</strong> : str</p>
<blockquote class="last">
<div><p>A filename to save to</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.set_t" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the value of the time of the world.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>t</strong> : Real</p>
<blockquote class="last">
<div><p>The time of the world</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.shape">
<code class="descname">shape</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Integer3<a class="headerlink" href="#ecell4.lattice.LatticeWorld.shape" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the triplet of sizes of column, row and layer.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.size">
<code class="descname">size</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.size" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the size of voxels.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.update_particle">
<code class="descname">update_particle</code><span class="sig-paren">(</span><em>pid</em>, <em>p</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.update_particle" title="Permalink to this definition">¶</a></dt>
<dd><p>Update a particle.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pid</strong> : ParticleID</p>
<blockquote>
<div><p>A particle id of the particle to update</p>
</div></blockquote>
<p><strong>p</strong> : Particle</p>
<blockquote>
<div><p>The information to update a particle</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>True if a new particle was created.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.update_voxel">
<code class="descname">update_voxel</code><span class="sig-paren">(</span><em>pid</em>, <em>v</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.lattice.LatticeWorld.update_voxel" title="Permalink to this definition">¶</a></dt>
<dd><p>Update a particle.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pid</strong> : ParticleID</p>
<blockquote>
<div><p>A particle id of the particle to update</p>
</div></blockquote>
<p><strong>v</strong> : Voxel</p>
<blockquote>
<div><p>The information to update</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>whether to succeed to update the particle</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.volume">
<code class="descname">volume</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the volume of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.voxel_radius">
<code class="descname">voxel_radius</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.voxel_radius" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the voxel radius.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.lattice.LatticeWorld.voxel_volume">
<code class="descname">voxel_volume</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.lattice.LatticeWorld.voxel_volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the volume of a voxel.</p>
</dd></dl>

</dd></dl>
