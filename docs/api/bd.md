  <span class="target" id="module-ecell4.bd"></span><dl class="class">
<dt id="ecell4.bd.BDFactory">
<em class="property">class </em><code class="descclassname">ecell4.bd.</code><code class="descname">BDFactory</code><a class="headerlink" href="#ecell4.bd.BDFactory" title="Permalink to this definition">¶</a></dt>
<dd><p>A factory class creating a BDWorld instance and a BDSimulator instance.</p>
<p>BDFactory(Integer3 matrix_sizes=None, GSLRandomNumberGenerator rng=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDFactory.create_simulator" title="ecell4.bd.BDFactory.create_simulator"><code class="xref py py-obj docutils literal"><span class="pre">create_simulator</span></code></a>((arg1[,&nbsp;arg2])</td>
<td>Return a <code class="docutils literal"><span class="pre">BDSimulator</span></code> instance.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDFactory.create_world" title="ecell4.bd.BDFactory.create_world"><code class="xref py py-obj docutils literal"><span class="pre">create_world</span></code></a>((arg1=None)&nbsp;-&gt;&nbsp;BDWorld)</td>
<td>Return a <code class="docutils literal"><span class="pre">BDWorld</span></code> instance.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.bd.BDFactory.create_simulator">
<code class="descname">create_simulator</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2=None</em><span class="sig-paren">)</span> &rarr; BDSimulator<a class="headerlink" href="#ecell4.bd.BDFactory.create_simulator" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a <code class="docutils literal"><span class="pre">BDSimulator</span></code> instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : BDWorld</p>
<blockquote>
<div><p>A world</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : Model</p>
<blockquote>
<div><p>A simulation model</p>
</div></blockquote>
<p><strong>arg2</strong> : BDWorld</p>
<blockquote>
<div><p>A world</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">BDSimulator:</p>
<blockquote class="last">
<div><p>The created simulator</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDFactory.create_world">
<code class="descname">create_world</code><span class="sig-paren">(</span><em>arg1=None</em><span class="sig-paren">)</span> &rarr; BDWorld<a class="headerlink" href="#ecell4.bd.BDFactory.create_world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a <code class="docutils literal"><span class="pre">BDWorld</span></code> instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : Real3</p>
<blockquote>
<div><p>The lengths of edges of a <code class="docutils literal"><span class="pre">BDWorld</span></code> created</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : str</p>
<blockquote>
<div><p>The path of a HDF5 file for <code class="docutils literal"><span class="pre">BDWorld</span></code></p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">BDWorld:</p>
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
<dt id="ecell4.bd.BDSimulator">
<em class="property">class </em><code class="descclassname">ecell4.bd.</code><code class="descname">BDSimulator</code><a class="headerlink" href="#ecell4.bd.BDSimulator" title="Permalink to this definition">¶</a></dt>
<dd><p>A class running the simulation with the bd algorithm.</p>
<p>BDSimulator(m, w)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.dt" title="ecell4.bd.BDSimulator.dt"><code class="xref py py-obj docutils literal"><span class="pre">dt</span></code></a></td>
<td>Return the step interval.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.initialize" title="ecell4.bd.BDSimulator.initialize"><code class="xref py py-obj docutils literal"><span class="pre">initialize</span></code></a></td>
<td>Initialize the simulator.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.last_reactions" title="ecell4.bd.BDSimulator.last_reactions"><code class="xref py py-obj docutils literal"><span class="pre">last_reactions</span></code></a>(()&nbsp;-&gt;&nbsp;[(ReactionRule,&nbsp;...)</td>
<td>Return reactions occuring at the last step.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.model" title="ecell4.bd.BDSimulator.model"><code class="xref py py-obj docutils literal"><span class="pre">model</span></code></a></td>
<td>Return the model bound.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.next_time" title="ecell4.bd.BDSimulator.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the scheduled time for the next step.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.num_steps" title="ecell4.bd.BDSimulator.num_steps"><code class="xref py py-obj docutils literal"><span class="pre">num_steps</span></code></a></td>
<td>Return the number of steps.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.run" title="ecell4.bd.BDSimulator.run"><code class="xref py py-obj docutils literal"><span class="pre">run</span></code></a>(duration,&nbsp;observers)</td>
<td>Run the simulation.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.set_dt" title="ecell4.bd.BDSimulator.set_dt"><code class="xref py py-obj docutils literal"><span class="pre">set_dt</span></code></a>(dt)</td>
<td>Set a step interval.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.set_t" title="ecell4.bd.BDSimulator.set_t"><code class="xref py py-obj docutils literal"><span class="pre">set_t</span></code></a>(t)</td>
<td>Set the current time.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.step" title="ecell4.bd.BDSimulator.step"><code class="xref py py-obj docutils literal"><span class="pre">step</span></code></a>((upto=None)&nbsp;-&gt;&nbsp;bool)</td>
<td>Step the simulation.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.t" title="ecell4.bd.BDSimulator.t"><code class="xref py py-obj docutils literal"><span class="pre">t</span></code></a></td>
<td>Return the time.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDSimulator.world" title="ecell4.bd.BDSimulator.world"><code class="xref py py-obj docutils literal"><span class="pre">world</span></code></a></td>
<td>Return the world bound.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.bd.BDSimulator.dt">
<code class="descname">dt</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDSimulator.dt" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the step interval.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDSimulator.initialize">
<code class="descname">initialize</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDSimulator.initialize" title="Permalink to this definition">¶</a></dt>
<dd><p>Initialize the simulator.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDSimulator.last_reactions">
<code class="descname">last_reactions</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [(ReactionRule, ReactionInfo)]<a class="headerlink" href="#ecell4.bd.BDSimulator.last_reactions" title="Permalink to this definition">¶</a></dt>
<dd><p>Return reactions occuring at the last step.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>the list of reaction rules and infos.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDSimulator.model">
<code class="descname">model</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDSimulator.model" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the model bound.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDSimulator.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDSimulator.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the scheduled time for the next step.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDSimulator.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDSimulator.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDSimulator.run">
<code class="descname">run</code><span class="sig-paren">(</span><em>duration</em>, <em>observers</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDSimulator.run" title="Permalink to this definition">¶</a></dt>
<dd><p>Run the simulation.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>duration</strong> : Real</p>
<blockquote>
<div><p>a duration for running a simulation.
A simulation is expected to be stopped at <code class="docutils literal"><span class="pre">t()</span> <span class="pre">+</span> <span class="pre">duration</span></code>.</p>
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
<dt id="ecell4.bd.BDSimulator.set_dt">
<code class="descname">set_dt</code><span class="sig-paren">(</span><em>dt</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDSimulator.set_dt" title="Permalink to this definition">¶</a></dt>
<dd><p>Set a step interval.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>dt</strong> : Real</p>
<blockquote class="last">
<div><p>a step interval</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDSimulator.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDSimulator.set_t" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the current time.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>t</strong> : Real</p>
<blockquote class="last">
<div><p>a current time.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDSimulator.step">
<code class="descname">step</code><span class="sig-paren">(</span><em>upto=None</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.bd.BDSimulator.step" title="Permalink to this definition">¶</a></dt>
<dd><p>Step the simulation.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>upto</strong> : Real, optional</p>
<blockquote>
<div><p>the time which to step the simulation up to</p>
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
<dt id="ecell4.bd.BDSimulator.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDSimulator.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDSimulator.world">
<code class="descname">world</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDSimulator.world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the world bound.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.bd.BDWorld">
<em class="property">class </em><code class="descclassname">ecell4.bd.</code><code class="descname">BDWorld</code><a class="headerlink" href="#ecell4.bd.BDWorld" title="Permalink to this definition">¶</a></dt>
<dd><p>A class containing the properties of the bd world.</p>
<p>BDWorld(edge_lengths=None, matrix_sizes=None, GSLRandomNumberGenerator rng=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.add_molecules" title="ecell4.bd.BDWorld.add_molecules"><code class="xref py py-obj docutils literal"><span class="pre">add_molecules</span></code></a>(sp,&nbsp;num[,&nbsp;shape])</td>
<td>Add some molecules.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.apply_boundary" title="ecell4.bd.BDWorld.apply_boundary"><code class="xref py py-obj docutils literal"><span class="pre">apply_boundary</span></code></a>((Real3&nbsp;pos)&nbsp;-&gt;&nbsp;Real3)</td>
<td>Return a position within the world by applying periodic boundaries to the given position.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.as_base" title="ecell4.bd.BDWorld.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Return self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.bind_to" title="ecell4.bd.BDWorld.bind_to"><code class="xref py py-obj docutils literal"><span class="pre">bind_to</span></code></a>(m)</td>
<td>Bind a model to the world</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.distance" title="ecell4.bd.BDWorld.distance"><code class="xref py py-obj docutils literal"><span class="pre">distance</span></code></a>((Real3&nbsp;pos1,&nbsp;Real3&nbsp;pos2)&nbsp;-&gt;&nbsp;Real)</td>
<td>Return the closest distance between the given positions.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.distance_sq" title="ecell4.bd.BDWorld.distance_sq"><code class="xref py py-obj docutils literal"><span class="pre">distance_sq</span></code></a>((Real3&nbsp;pos1,&nbsp;Real3&nbsp;pos2)&nbsp;-&gt;&nbsp;Real)</td>
<td>Return a square of the closest distance between the given positions.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.edge_lengths" title="ecell4.bd.BDWorld.edge_lengths"><code class="xref py py-obj docutils literal"><span class="pre">edge_lengths</span></code></a>(()&nbsp;-&gt;&nbsp;Real3)</td>
<td>Return the edge lengths of the world.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.get_particle" title="ecell4.bd.BDWorld.get_particle"><code class="xref py py-obj docutils literal"><span class="pre">get_particle</span></code></a>(pid)&nbsp;-&gt;&nbsp;(ParticleID,&nbsp;Particle)</td>
<td>Return the particle associated a given ParticleID.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.has_particle" title="ecell4.bd.BDWorld.has_particle"><code class="xref py py-obj docutils literal"><span class="pre">has_particle</span></code></a>((pid)&nbsp;-&gt;&nbsp;bool)</td>
<td>Check if a particle associated with a given particle id exists.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.list_particles" title="ecell4.bd.BDWorld.list_particles"><code class="xref py py-obj docutils literal"><span class="pre">list_particles</span></code></a>((sp)&nbsp;-&gt;&nbsp;[(ParticleID,&nbsp;Particle)])</td>
<td>Return the list of particles.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.list_particles_exact" title="ecell4.bd.BDWorld.list_particles_exact"><code class="xref py py-obj docutils literal"><span class="pre">list_particles_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;[(ParticleID,&nbsp;...)</td>
<td>Return the list of particles of a given species.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.list_particles_within_radius" title="ecell4.bd.BDWorld.list_particles_within_radius"><code class="xref py py-obj docutils literal"><span class="pre">list_particles_within_radius</span></code></a>((pos,&nbsp;radius[,&nbsp;...])</td>
<td>Returns a list of pairs of a particle and distance within the given sphere.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.load" title="ecell4.bd.BDWorld.load"><code class="xref py py-obj docutils literal"><span class="pre">load</span></code></a>(filename)</td>
<td>Load the world from a file.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.new_particle" title="ecell4.bd.BDWorld.new_particle"><code class="xref py py-obj docutils literal"><span class="pre">new_particle</span></code></a>(arg1[,&nbsp;arg2])</td>
<td>Create a new particle.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.num_molecules" title="ecell4.bd.BDWorld.num_molecules"><code class="xref py py-obj docutils literal"><span class="pre">num_molecules</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of molecules.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.num_molecules_exact" title="ecell4.bd.BDWorld.num_molecules_exact"><code class="xref py py-obj docutils literal"><span class="pre">num_molecules_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of molecules of a given species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.num_particles" title="ecell4.bd.BDWorld.num_particles"><code class="xref py py-obj docutils literal"><span class="pre">num_particles</span></code></a>((sp=None)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of particles.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.num_particles_exact" title="ecell4.bd.BDWorld.num_particles_exact"><code class="xref py py-obj docutils literal"><span class="pre">num_particles_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of particles of a given species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.periodic_transpose" title="ecell4.bd.BDWorld.periodic_transpose"><code class="xref py py-obj docutils literal"><span class="pre">periodic_transpose</span></code></a>((Real3&nbsp;pos1,&nbsp;...)</td>
<td>Return a closest image of pos1 relative to the given position (pos2).</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.remove_molecules" title="ecell4.bd.BDWorld.remove_molecules"><code class="xref py py-obj docutils literal"><span class="pre">remove_molecules</span></code></a>(sp,&nbsp;num)</td>
<td>Remove the molecules.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.remove_particle" title="ecell4.bd.BDWorld.remove_particle"><code class="xref py py-obj docutils literal"><span class="pre">remove_particle</span></code></a>(pid)</td>
<td>Remove the particle associated with a given ParticleID.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.rng" title="ecell4.bd.BDWorld.rng"><code class="xref py py-obj docutils literal"><span class="pre">rng</span></code></a></td>
<td>Return a random number generator object.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.save" title="ecell4.bd.BDWorld.save"><code class="xref py py-obj docutils literal"><span class="pre">save</span></code></a>(filename)</td>
<td>Save the world to a file.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.set_t" title="ecell4.bd.BDWorld.set_t"><code class="xref py py-obj docutils literal"><span class="pre">set_t</span></code></a>(t)</td>
<td>Set the value of the time of the world.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.t" title="ecell4.bd.BDWorld.t"><code class="xref py py-obj docutils literal"><span class="pre">t</span></code></a></td>
<td>Return the time of the world.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.bd.BDWorld.update_particle" title="ecell4.bd.BDWorld.update_particle"><code class="xref py py-obj docutils literal"><span class="pre">update_particle</span></code></a>((pid,&nbsp;p)&nbsp;-&gt;&nbsp;bool)</td>
<td>Update a particle.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.bd.BDWorld.volume" title="ecell4.bd.BDWorld.volume"><code class="xref py py-obj docutils literal"><span class="pre">volume</span></code></a></td>
<td>Return the volume of the world.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.bd.BDWorld.add_molecules">
<code class="descname">add_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em>, <em>shape=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.add_molecules" title="Permalink to this definition">¶</a></dt>
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
<p><strong>shape</strong> : Shape, optional</p>
<blockquote class="last">
<div><p>a shape to add molecules on [not supported yet]</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.apply_boundary">
<code class="descname">apply_boundary</code><span class="sig-paren">(</span><em>Real3 pos</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.bd.BDWorld.apply_boundary" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a position within the world by applying periodic boundaries
to the given position.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Return self as a base class. Only for developmental use.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.bind_to">
<code class="descname">bind_to</code><span class="sig-paren">(</span><em>m</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.bind_to" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.distance">
<code class="descname">distance</code><span class="sig-paren">(</span><em>Real3 pos1</em>, <em>Real3 pos2</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.bd.BDWorld.distance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the closest distance between the given positions.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.distance_sq">
<code class="descname">distance_sq</code><span class="sig-paren">(</span><em>Real3 pos1</em>, <em>Real3 pos2</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.bd.BDWorld.distance_sq" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a square of the closest distance between the given positions.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.edge_lengths">
<code class="descname">edge_lengths</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.bd.BDWorld.edge_lengths" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the edge lengths of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.get_particle">
<code class="descname">get_particle</code><span class="sig-paren">(</span><em>pid) -&gt; (ParticleID</em>, <em>Particle</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.get_particle" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the particle associated a given ParticleID.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pid</strong> : ParticleID</p>
<blockquote>
<div><p>An id of the particle you want</p>
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
<dt id="ecell4.bd.BDWorld.has_particle">
<code class="descname">has_particle</code><span class="sig-paren">(</span><em>pid</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.bd.BDWorld.has_particle" title="Permalink to this definition">¶</a></dt>
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
<div><p>If a particle exists, return True. Otherwise return False</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.list_particles">
<code class="descname">list_particles</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.bd.BDWorld.list_particles" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.list_particles_exact">
<code class="descname">list_particles_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.bd.BDWorld.list_particles_exact" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.list_particles_within_radius">
<code class="descname">list_particles_within_radius</code><span class="sig-paren">(</span><em>pos</em>, <em>radius</em>, <em>ignore1=None</em>, <em>ignore2=None</em><span class="sig-paren">)</span> &rarr; [((ParticleID, Particle), Real)]<a class="headerlink" href="#ecell4.bd.BDWorld.list_particles_within_radius" title="Permalink to this definition">¶</a></dt>
<dd><p>Returns a list of pairs of a particle and distance within the given sphere.
The region is specified with a center position and radius.
ignore1 and ignore2 will be removed from the list.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A center position.</p>
</div></blockquote>
<p><strong>radius</strong> : Real</p>
<blockquote>
<div><p>A radius.</p>
</div></blockquote>
<p><strong>ignore1</strong> : ParticleID, optional</p>
<blockquote>
<div><p>An id ignored.</p>
</div></blockquote>
<p><strong>ignore2</strong> : ParticleID, optional</p>
<blockquote>
<div><p>An id ignored.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of pairs of a particle and its distance from the center position.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.load">
<code class="descname">load</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.load" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.new_particle">
<code class="descname">new_particle</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2=None) -&gt; (ParticleID</em>, <em>Particle</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.new_particle" title="Permalink to this definition">¶</a></dt>
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
<div><p>A position to place a particle</p>
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
<dt id="ecell4.bd.BDWorld.num_molecules">
<code class="descname">num_molecules</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.bd.BDWorld.num_molecules" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.num_molecules_exact">
<code class="descname">num_molecules_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.bd.BDWorld.num_molecules_exact" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.num_particles">
<code class="descname">num_particles</code><span class="sig-paren">(</span><em>sp=None</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.bd.BDWorld.num_particles" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.num_particles_exact">
<code class="descname">num_particles_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.bd.BDWorld.num_particles_exact" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.periodic_transpose">
<code class="descname">periodic_transpose</code><span class="sig-paren">(</span><em>Real3 pos1</em>, <em>Real3 pos2</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.bd.BDWorld.periodic_transpose" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a closest image of pos1 relative to the given position (pos2).</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.remove_molecules">
<code class="descname">remove_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.remove_molecules" title="Permalink to this definition">¶</a></dt>
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
<blockquote class="last">
<div><p>a number of molecules to be removed</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.remove_particle">
<code class="descname">remove_particle</code><span class="sig-paren">(</span><em>pid</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.remove_particle" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove the particle associated with a given ParticleID.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pid</strong> : ParticleID</p>
<blockquote class="last">
<div><p>An id of particle to remove</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.rng">
<code class="descname">rng</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.rng" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a random number generator object.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.save">
<code class="descname">save</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.save" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.set_t" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.bd.BDWorld.update_particle">
<code class="descname">update_particle</code><span class="sig-paren">(</span><em>pid</em>, <em>p</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.bd.BDWorld.update_particle" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.bd.BDWorld.volume">
<code class="descname">volume</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.bd.BDWorld.volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the volume of the world.</p>
</dd></dl>

</dd></dl>
