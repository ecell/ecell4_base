  <span class="target" id="module-ecell4.gillespie"></span><dl class="class">
<dt id="ecell4.gillespie.GillespieFactory">
<em class="property">class </em><code class="descclassname">ecell4.gillespie.</code><code class="descname">GillespieFactory</code><a class="headerlink" href="#ecell4.gillespie.GillespieFactory" title="Permalink to this definition">¶</a></dt>
<dd><p>A factory class creating a GillespieWorld instance and a GillespieSimulator instance.</p>
<p>GillespieFactory(GSLRandomNumberGenerator rng=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieFactory.create_simulator" title="ecell4.gillespie.GillespieFactory.create_simulator"><code class="xref py py-obj docutils literal"><span class="pre">create_simulator</span></code></a>((arg1,&nbsp;...)</td>
<td>Return a GillespieSimulator instance.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieFactory.create_world" title="ecell4.gillespie.GillespieFactory.create_world"><code class="xref py py-obj docutils literal"><span class="pre">create_world</span></code></a>((arg1=None)&nbsp;-&gt;&nbsp;GillespieWorld)</td>
<td>Return a GillespieWorld instance.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.gillespie.GillespieFactory.create_simulator">
<code class="descname">create_simulator</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2</em><span class="sig-paren">)</span> &rarr; GillespieSimulator<a class="headerlink" href="#ecell4.gillespie.GillespieFactory.create_simulator" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a GillespieSimulator instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : GillespieWorld</p>
<blockquote>
<div><p>a world</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : Model</p>
<blockquote>
<div><p>a simulation model</p>
</div></blockquote>
<p><strong>arg2</strong> : GillespieWorld</p>
<blockquote>
<div><p>a world</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">GillespieSimulator:</p>
<blockquote class="last">
<div><p>the created simulator</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieFactory.create_world">
<code class="descname">create_world</code><span class="sig-paren">(</span><em>arg1=None</em><span class="sig-paren">)</span> &rarr; GillespieWorld<a class="headerlink" href="#ecell4.gillespie.GillespieFactory.create_world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a GillespieWorld instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : Real3</p>
<blockquote>
<div><p>The lengths of edges of a GillespieWorld created</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : str</p>
<blockquote>
<div><p>The path of a HDF5 file for GillespieWorld</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">GillespieWorld:</p>
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
<dt id="ecell4.gillespie.GillespieSimulator">
<em class="property">class </em><code class="descclassname">ecell4.gillespie.</code><code class="descname">GillespieSimulator</code><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator" title="Permalink to this definition">¶</a></dt>
<dd><p>A class running the simulation with the gillespie algorithm.</p>
<p>GillespieSimulator(m, w)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.dt" title="ecell4.gillespie.GillespieSimulator.dt"><code class="xref py py-obj docutils literal"><span class="pre">dt</span></code></a></td>
<td>Return the step interval.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.initialize" title="ecell4.gillespie.GillespieSimulator.initialize"><code class="xref py py-obj docutils literal"><span class="pre">initialize</span></code></a></td>
<td>Initialize the simulator.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.last_reactions" title="ecell4.gillespie.GillespieSimulator.last_reactions"><code class="xref py py-obj docutils literal"><span class="pre">last_reactions</span></code></a>(()&nbsp;-&gt;&nbsp;[(ReactionRule,&nbsp;...)</td>
<td>Return reactions occuring at the last step.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.model" title="ecell4.gillespie.GillespieSimulator.model"><code class="xref py py-obj docutils literal"><span class="pre">model</span></code></a></td>
<td>Return the model bound.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.next_time" title="ecell4.gillespie.GillespieSimulator.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the scheduled time for the next step.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.num_steps" title="ecell4.gillespie.GillespieSimulator.num_steps"><code class="xref py py-obj docutils literal"><span class="pre">num_steps</span></code></a></td>
<td>Return the number of steps.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.run" title="ecell4.gillespie.GillespieSimulator.run"><code class="xref py py-obj docutils literal"><span class="pre">run</span></code></a>(duration,&nbsp;observers)</td>
<td>Run the simulation.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.set_dt" title="ecell4.gillespie.GillespieSimulator.set_dt"><code class="xref py py-obj docutils literal"><span class="pre">set_dt</span></code></a>(dt)</td>
<td>Set a step interval.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.set_t" title="ecell4.gillespie.GillespieSimulator.set_t"><code class="xref py py-obj docutils literal"><span class="pre">set_t</span></code></a>(t)</td>
<td>Set the current time.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.step" title="ecell4.gillespie.GillespieSimulator.step"><code class="xref py py-obj docutils literal"><span class="pre">step</span></code></a>((upto=None)&nbsp;-&gt;&nbsp;bool)</td>
<td>Step the simulation.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.t" title="ecell4.gillespie.GillespieSimulator.t"><code class="xref py py-obj docutils literal"><span class="pre">t</span></code></a></td>
<td>Return the time.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieSimulator.world" title="ecell4.gillespie.GillespieSimulator.world"><code class="xref py py-obj docutils literal"><span class="pre">world</span></code></a></td>
<td>Return the world bound.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.gillespie.GillespieSimulator.dt">
<code class="descname">dt</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.dt" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the step interval.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieSimulator.initialize">
<code class="descname">initialize</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.initialize" title="Permalink to this definition">¶</a></dt>
<dd><p>Initialize the simulator.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieSimulator.last_reactions">
<code class="descname">last_reactions</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [(ReactionRule, ReactionInfo)]<a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.last_reactions" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.gillespie.GillespieSimulator.model">
<code class="descname">model</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.model" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the model bound.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieSimulator.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the scheduled time for the next step.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieSimulator.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieSimulator.run">
<code class="descname">run</code><span class="sig-paren">(</span><em>duration</em>, <em>observers</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.run" title="Permalink to this definition">¶</a></dt>
<dd><p>Run the simulation.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>duration</strong> : Real</p>
<blockquote>
<div><dl class="docutils">
<dt>a duration for running a simulation.</dt>
<dd><p class="first last">A simulation is expected to be stopped at t() + duration.</p>
</dd>
</dl>
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
<dt id="ecell4.gillespie.GillespieSimulator.set_dt">
<code class="descname">set_dt</code><span class="sig-paren">(</span><em>dt</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.set_dt" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.gillespie.GillespieSimulator.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.set_t" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.gillespie.GillespieSimulator.step">
<code class="descname">step</code><span class="sig-paren">(</span><em>upto=None</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.step" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.gillespie.GillespieSimulator.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieSimulator.world">
<code class="descname">world</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieSimulator.world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the world bound.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.gillespie.GillespieWorld">
<em class="property">class </em><code class="descclassname">ecell4.gillespie.</code><code class="descname">GillespieWorld</code><a class="headerlink" href="#ecell4.gillespie.GillespieWorld" title="Permalink to this definition">¶</a></dt>
<dd><p>A class containing the properties of the gillespie world.</p>
<p>GillespieWorld(edge_lengths=None, GSLRandomNumberGenerator rng=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.add_molecules" title="ecell4.gillespie.GillespieWorld.add_molecules"><code class="xref py py-obj docutils literal"><span class="pre">add_molecules</span></code></a>(sp,&nbsp;num[,&nbsp;shape])</td>
<td>Add some molecules.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.as_base" title="ecell4.gillespie.GillespieWorld.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Return self as a base class.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.bind_to" title="ecell4.gillespie.GillespieWorld.bind_to"><code class="xref py py-obj docutils literal"><span class="pre">bind_to</span></code></a>(m)</td>
<td>Bind a model to the world</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.edge_lengths" title="ecell4.gillespie.GillespieWorld.edge_lengths"><code class="xref py py-obj docutils literal"><span class="pre">edge_lengths</span></code></a>(()&nbsp;-&gt;&nbsp;Real3)</td>
<td>Return the edge lengths of the world.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.list_particles" title="ecell4.gillespie.GillespieWorld.list_particles"><code class="xref py py-obj docutils literal"><span class="pre">list_particles</span></code></a>((sp)&nbsp;-&gt;&nbsp;[(ParticleID,&nbsp;Particle)])</td>
<td>Return the list of particles.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.list_particles_exact" title="ecell4.gillespie.GillespieWorld.list_particles_exact"><code class="xref py py-obj docutils literal"><span class="pre">list_particles_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;[(ParticleID,&nbsp;...)</td>
<td>Return the list of particles of a given species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.list_species" title="ecell4.gillespie.GillespieWorld.list_species"><code class="xref py py-obj docutils literal"><span class="pre">list_species</span></code></a></td>
<td>Return a list of species.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.load" title="ecell4.gillespie.GillespieWorld.load"><code class="xref py py-obj docutils literal"><span class="pre">load</span></code></a>(filename)</td>
<td>Load self from a HDF5 file.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.num_molecules" title="ecell4.gillespie.GillespieWorld.num_molecules"><code class="xref py py-obj docutils literal"><span class="pre">num_molecules</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of molecules.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.num_molecules_exact" title="ecell4.gillespie.GillespieWorld.num_molecules_exact"><code class="xref py py-obj docutils literal"><span class="pre">num_molecules_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of molecules of a given species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.remove_molecules" title="ecell4.gillespie.GillespieWorld.remove_molecules"><code class="xref py py-obj docutils literal"><span class="pre">remove_molecules</span></code></a>(sp,&nbsp;num)</td>
<td>Remove the molecules.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.rng" title="ecell4.gillespie.GillespieWorld.rng"><code class="xref py py-obj docutils literal"><span class="pre">rng</span></code></a></td>
<td>Return a random number generator object.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.save" title="ecell4.gillespie.GillespieWorld.save"><code class="xref py py-obj docutils literal"><span class="pre">save</span></code></a>(filename)</td>
<td>Save self to a HDF5 file.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.set_t" title="ecell4.gillespie.GillespieWorld.set_t"><code class="xref py py-obj docutils literal"><span class="pre">set_t</span></code></a>(t)</td>
<td>Set the value of the time of the world.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.t" title="ecell4.gillespie.GillespieWorld.t"><code class="xref py py-obj docutils literal"><span class="pre">t</span></code></a></td>
<td>Return the time of the world.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.gillespie.GillespieWorld.volume" title="ecell4.gillespie.GillespieWorld.volume"><code class="xref py py-obj docutils literal"><span class="pre">volume</span></code></a></td>
<td>Return the volume of the world.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.add_molecules">
<code class="descname">add_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em>, <em>shape=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.add_molecules" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.gillespie.GillespieWorld.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Return self as a base class. Only for developmental use.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.bind_to">
<code class="descname">bind_to</code><span class="sig-paren">(</span><em>m</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.bind_to" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.gillespie.GillespieWorld.edge_lengths">
<code class="descname">edge_lengths</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.gillespie.GillespieWorld.edge_lengths" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the edge lengths of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.list_particles">
<code class="descname">list_particles</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.gillespie.GillespieWorld.list_particles" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the list of particles.
A position of each particle is randomly generated.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species, optional</p>
<blockquote>
<div><dl class="docutils">
<dt>the species of particles to list up</dt>
<dd><p class="first last">If no species is given, return the whole list of particles.</p>
</dd>
</dl>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>the list of particles (of the given species)</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.list_particles_exact">
<code class="descname">list_particles_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.gillespie.GillespieWorld.list_particles_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the list of particles of a given species.
A position of each particle is randomly generated.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>the species of particles to list up</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>the list of particles of a given species</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.list_species">
<code class="descname">list_species</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.list_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of species.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.load">
<code class="descname">load</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.load" title="Permalink to this definition">¶</a></dt>
<dd><p>Load self from a HDF5 file.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>filename</strong> : str</p>
<blockquote class="last">
<div><p>a filename</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.num_molecules">
<code class="descname">num_molecules</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.gillespie.GillespieWorld.num_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of molecules.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>a species whose molecules you count</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>the number of molecules (of a given species)</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.num_molecules_exact">
<code class="descname">num_molecules_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.gillespie.GillespieWorld.num_molecules_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of molecules of a given species.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>a species whose molecules you count</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>the number of molecules of a given species</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.remove_molecules">
<code class="descname">remove_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.remove_molecules" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.gillespie.GillespieWorld.rng">
<code class="descname">rng</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.rng" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a random number generator object.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.save">
<code class="descname">save</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.save" title="Permalink to this definition">¶</a></dt>
<dd><p>Save self to a HDF5 file.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>filename</strong> : str</p>
<blockquote class="last">
<div><p>a filename</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.set_t" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the value of the time of the world.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>t</strong> : Real</p>
<blockquote class="last">
<div><p>the time of the world</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.gillespie.GillespieWorld.volume">
<code class="descname">volume</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.gillespie.GillespieWorld.volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the volume of the world.</p>
</dd></dl>

</dd></dl>
