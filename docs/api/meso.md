  <span class="target" id="module-ecell4.meso"></span><dl class="class">
<dt id="ecell4.meso.MesoscopicFactory">
<em class="property">class </em><code class="descclassname">ecell4.meso.</code><code class="descname">MesoscopicFactory</code><a class="headerlink" href="#ecell4.meso.MesoscopicFactory" title="Permalink to this definition">¶</a></dt>
<dd><p>A factory class creating a MesoscopicWorld instance and a MesoscopicSimulator instance.</p>
<p>MesoscopicFactory(matrix_sizes=None, GSLRandomNumberGenerator rng=None)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.meso.MesoscopicFactory.create_simulator">
<code class="descname">create_simulator</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2</em><span class="sig-paren">)</span> &rarr; MesoscopicSimulator<a class="headerlink" href="#ecell4.meso.MesoscopicFactory.create_simulator" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a MesoscopicSimulator instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : MesoscopicWorld</p>
<blockquote>
<div><p>a world</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : Model</p>
<blockquote>
<div><p>a simulation model</p>
</div></blockquote>
<p><strong>arg2</strong> : MesoscopicWorld</p>
<blockquote>
<div><p>a world</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">MesoscopicSimulator:</p>
<blockquote class="last">
<div><p>the created simulator</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicFactory.create_world">
<code class="descname">create_world</code><span class="sig-paren">(</span><em>arg1=None</em><span class="sig-paren">)</span> &rarr; MesoscopicWorld<a class="headerlink" href="#ecell4.meso.MesoscopicFactory.create_world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a MesoscopicWorld instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : Real3</p>
<blockquote>
<div><p>The lengths of edges of a MesoscopicWorld created</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : str</p>
<blockquote>
<div><p>The path of a HDF5 file for MesoscopicWorld</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">MesoscopicWorld:</p>
<blockquote class="last">
<div><p>the created world</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.meso.MesoscopicSimulator">
<em class="property">class </em><code class="descclassname">ecell4.meso.</code><code class="descname">MesoscopicSimulator</code><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator" title="Permalink to this definition">¶</a></dt>
<dd><p>A class running the simulation with the meso algorithm.</p>
<p>MesoscopicSimulator(m, w)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.dt">
<code class="descname">dt</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.dt" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the step interval.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.initialize">
<code class="descname">initialize</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.initialize" title="Permalink to this definition">¶</a></dt>
<dd><p>Initialize the simulator.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.last_reactions">
<code class="descname">last_reactions</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [(ReactionRule, ReactionInfo)]<a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.last_reactions" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.meso.MesoscopicSimulator.model">
<code class="descname">model</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.model" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the model bound.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the scheduled time for the next step.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.run">
<code class="descname">run</code><span class="sig-paren">(</span><em>duration</em>, <em>observers</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.run" title="Permalink to this definition">¶</a></dt>
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
<div><p>Observers</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.set_dt">
<code class="descname">set_dt</code><span class="sig-paren">(</span><em>dt</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.set_dt" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.meso.MesoscopicSimulator.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.set_t" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.meso.MesoscopicSimulator.step">
<code class="descname">step</code><span class="sig-paren">(</span><em>upto=None</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.step" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.meso.MesoscopicSimulator.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.world">
<code class="descname">world</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the world bound.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.meso.MesoscopicWorld">
<em class="property">class </em><code class="descclassname">ecell4.meso.</code><code class="descname">MesoscopicWorld</code><a class="headerlink" href="#ecell4.meso.MesoscopicWorld" title="Permalink to this definition">¶</a></dt>
<dd><p>A class containing the properties of the meso world.</p>
<p>MesoscopicWorld(edge_lengths=None, matrix_sizes=None, GSLRandomNumberGenerator rng=None)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.add_molecules">
<code class="descname">add_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em>, <em>c=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.add_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Add some molecules.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>a species of molecules to add</p>
</div></blockquote>
<p><strong>num</strong> : Integer</p>
<blockquote>
<div><p>the number of molecules to add</p>
</div></blockquote>
<p><strong>c</strong> : Integer or Shape, optional</p>
<blockquote class="last">
<div><p>a coordinate or shape to add molecules on</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.add_structure">
<code class="descname">add_structure</code><span class="sig-paren">(</span><em>sp</em>, <em>shape</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.add_structure" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a structure.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>a species suggesting the shape.</p>
</div></blockquote>
<p><strong>shape</strong> : Shape</p>
<blockquote class="last">
<div><p>a shape of the structure.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Return self as a base class. Only for developmental use.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.bind_to">
<code class="descname">bind_to</code><span class="sig-paren">(</span><em>m</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.bind_to" title="Permalink to this definition">¶</a></dt>
<dd><p>Bind a model to the world</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>m</strong> : Model</p>
<blockquote class="last">
<div><p>a model to bind</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.check_structure">
<code class="descname">check_structure</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.check_structure" title="Permalink to this definition">¶</a></dt>
<dd><p>on_structure(sp, g) -&gt; bool</p>
<p>Check if the given subvolume is belonging to the structure.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species for the target structure.</p>
</div></blockquote>
<p><strong>g</strong> : Integer3</p>
<blockquote>
<div><p>a global coordinate pointing a subvolume</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>True if the subvolume is belonging to the structure.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.coord2global">
<code class="descname">coord2global</code><span class="sig-paren">(</span><em>coord</em><span class="sig-paren">)</span> &rarr; Integer3<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.coord2global" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a coordinate to a global coordinate.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.edge_lengths">
<code class="descname">edge_lengths</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.edge_lengths" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the edge lengths of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.get_occupancy">
<code class="descname">get_occupancy</code><span class="sig-paren">(</span><em>sp</em>, <em>g</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.get_occupancy" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the occupancy of the structure in the subvolume.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species for the target structure.</p>
</div></blockquote>
<p><strong>g</strong> : Integer3</p>
<blockquote>
<div><p>a global coordinate pointing a subvolume</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real:</p>
<blockquote class="last">
<div><p>The occupancy of the structure.
As a default, return 1 if the subvolume overlaps with the structure,
and return 0 otherwise.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.get_volume">
<code class="descname">get_volume</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.get_volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a volume of the given structure.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species for the target structure.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real:</p>
<blockquote class="last">
<div><p>A total volume of subvolumes belonging to the structure.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.global2coord">
<code class="descname">global2coord</code><span class="sig-paren">(</span><em>g</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.global2coord" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a global coordinate to a coordinate.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>g</strong> : Integer3</p>
<blockquote>
<div><p>a global coordinate</p>
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
<dt id="ecell4.meso.MesoscopicWorld.has_structure">
<code class="descname">has_structure</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.has_structure" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given structure is in the space or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species for the target structure.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>True if the given structure is in self.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.list_coordinates">
<code class="descname">list_coordinates</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [Integer]<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.list_coordinates" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of coordinates of molecules belonging to the given species.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species of molecules.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of coordinates.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.list_coordinates_exact">
<code class="descname">list_coordinates_exact</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.list_coordinates_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>list_coordinates(sp) -&gt; [Integer]</p>
<p>Return a list of coordinates of molecules belonging to the given species.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species of molecules.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of coordinates.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.list_particles">
<code class="descname">list_particles</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.list_particles" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.meso.MesoscopicWorld.list_particles_exact">
<code class="descname">list_particles_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.list_particles_exact" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.meso.MesoscopicWorld.list_species">
<code class="descname">list_species</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.list_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of species.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.load">
<code class="descname">load</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.load" title="Permalink to this definition">¶</a></dt>
<dd><p>Load the world from a file.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>filename</strong> : str</p>
<blockquote class="last">
<div><p>a filename to load from</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.matrix_sizes">
<code class="descname">matrix_sizes</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Integer3<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.matrix_sizes" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of subvolumes along axes.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.num_molecules">
<code class="descname">num_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>c=None</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.num_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of molecules within the suggested subvolume.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species whose molecules you count</p>
</div></blockquote>
<p><strong>c</strong> : Integer, optional</p>
<blockquote>
<div><p>A coordinate.</p>
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
<dt id="ecell4.meso.MesoscopicWorld.num_molecules_exact">
<code class="descname">num_molecules_exact</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.num_molecules_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>num_particles_exact(sp, c=None) -&gt; Integer</p>
<p>Return the number of molecules within the suggested subvolume.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>The species of molecules to count</p>
</div></blockquote>
<p><strong>c</strong> : Integer, optional</p>
<blockquote>
<div><p>A coordinate.</p>
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
<dt id="ecell4.meso.MesoscopicWorld.num_subvolumes">
<code class="descname">num_subvolumes</code><span class="sig-paren">(</span><em>sp=None</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.num_subvolumes" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of subvolumes.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species, optional</p>
<blockquote>
<div><p>A species specifying a structure.
When no species is given, return the total number of subvolumes.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of subvolumes.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.on_structure">
<code class="descname">on_structure</code><span class="sig-paren">(</span><em>sp</em>, <em>g</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.on_structure" title="Permalink to this definition">¶</a></dt>
<dd><p>Check if the given species would be on the proper structure at the coordinate.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>a species scheduled to be placed</p>
</div></blockquote>
<p><strong>g</strong> : Integer3</p>
<blockquote>
<div><p>a global coordinate pointing a subvolume</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>True if it is on the proper structure.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.position2coordinate">
<code class="descname">position2coordinate</code><span class="sig-paren">(</span><em>pos</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.position2coordinate" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.meso.MesoscopicWorld.position2global">
<code class="descname">position2global</code><span class="sig-paren">(</span><em>pos</em><span class="sig-paren">)</span> &rarr; Integer3<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.position2global" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.meso.MesoscopicWorld.remove_molecules">
<code class="descname">remove_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em>, <em>c=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.remove_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove the molecules.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>a species whose molecules to remove</p>
</div></blockquote>
<p><strong>num</strong> : Integer</p>
<blockquote>
<div><p>a number of molecules to be removed</p>
</div></blockquote>
<p><strong>c</strong> : Integer, optional</p>
<blockquote class="last">
<div><p>A coordinate.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.rng">
<code class="descname">rng</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.rng" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a random number generator object.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.save">
<code class="descname">save</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.save" title="Permalink to this definition">¶</a></dt>
<dd><p>Save the world to a file.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>filename</strong> : str</p>
<blockquote class="last">
<div><p>a filename to save to</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.set_t" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.meso.MesoscopicWorld.subvolume">
<code class="descname">subvolume</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.subvolume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the subvolume of each cell.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.subvolume_edge_lengths">
<code class="descname">subvolume_edge_lengths</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.subvolume_edge_lengths" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the edge lengths of a subvolume.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.volume">
<code class="descname">volume</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the volume of the world.</p>
</dd></dl>

</dd></dl>
