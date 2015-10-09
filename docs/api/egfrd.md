  <span class="target" id="module-ecell4.egfrd"></span><dl class="class">
<dt id="ecell4.egfrd.BDFactory">
<em class="property">class </em><code class="descclassname">ecell4.egfrd.</code><code class="descname">BDFactory</code><a class="headerlink" href="#ecell4.egfrd.BDFactory" title="Permalink to this definition">¶</a></dt>
<dd><p>A factory class creating a BDWorld instance and a BDSimulator instance.</p>
<dl class="docutils">
<dt>BDFactory(matrix_sizes=None, rng=None, dissociation_retry_moves=None,</dt>
<dd>bd_dt_factor=None)</dd>
</dl>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.egfrd.BDFactory.create_simulator">
<code class="descname">create_simulator</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2</em><span class="sig-paren">)</span> &rarr; BDSimulator<a class="headerlink" href="#ecell4.egfrd.BDFactory.create_simulator" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a BDSimulator instance.</p>
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
<dt id="ecell4.egfrd.BDFactory.create_world">
<code class="descname">create_world</code><span class="sig-paren">(</span><em>arg1=None</em><span class="sig-paren">)</span> &rarr; EGFRDWorld<a class="headerlink" href="#ecell4.egfrd.BDFactory.create_world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a EGFRDWorld instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : Real3</p>
<blockquote>
<div><p>The lengths of edges of a EGFRDWorld created</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : str</p>
<blockquote>
<div><p>The path of a HDF5 file for EGFRDWorld</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">EGFRDWorld:</p>
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
<dt id="ecell4.egfrd.BDSimulator">
<em class="property">class </em><code class="descclassname">ecell4.egfrd.</code><code class="descname">BDSimulator</code><a class="headerlink" href="#ecell4.egfrd.BDSimulator" title="Permalink to this definition">¶</a></dt>
<dd><p>A class running the simulation with the bd algorithm.</p>
<p>BDSimulator(m, w)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.egfrd.BDSimulator.dt">
<code class="descname">dt</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.BDSimulator.dt" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the step interval.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.BDSimulator.initialize">
<code class="descname">initialize</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.BDSimulator.initialize" title="Permalink to this definition">¶</a></dt>
<dd><p>Initialize the simulator.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.BDSimulator.last_reactions">
<code class="descname">last_reactions</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [(ReactionRule, ReactionInfo)]<a class="headerlink" href="#ecell4.egfrd.BDSimulator.last_reactions" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.BDSimulator.model">
<code class="descname">model</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.BDSimulator.model" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the model bound.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.BDSimulator.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.BDSimulator.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the scheduled time for the next step.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.BDSimulator.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.BDSimulator.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.BDSimulator.run">
<code class="descname">run</code><span class="sig-paren">(</span><em>duration</em>, <em>observers</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.BDSimulator.run" title="Permalink to this definition">¶</a></dt>
<dd><p>Run the simulation.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>duration</strong> : Real</p>
<blockquote>
<div><p>a duration for running a simulation.
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
<dt id="ecell4.egfrd.BDSimulator.set_dt">
<code class="descname">set_dt</code><span class="sig-paren">(</span><em>dt</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.BDSimulator.set_dt" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.BDSimulator.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.BDSimulator.set_t" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.BDSimulator.step">
<code class="descname">step</code><span class="sig-paren">(</span><em>upto=None</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.egfrd.BDSimulator.step" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.BDSimulator.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.BDSimulator.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.BDSimulator.world">
<code class="descname">world</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.BDSimulator.world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the world bound.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.egfrd.EGFRDFactory">
<em class="property">class </em><code class="descclassname">ecell4.egfrd.</code><code class="descname">EGFRDFactory</code><a class="headerlink" href="#ecell4.egfrd.EGFRDFactory" title="Permalink to this definition">¶</a></dt>
<dd><p>A factory class creating a BDWorld instance and a BDSimulator instance.</p>
<dl class="docutils">
<dt>EGFRDFactory(matrix_sizes=None, rng=None, dissociation_retry_moves,</dt>
<dd>bd_dt_factor, user_max_shell_size)</dd>
</dl>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.egfrd.EGFRDFactory.create_simulator">
<code class="descname">create_simulator</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2</em><span class="sig-paren">)</span> &rarr; EGFRDSimulator<a class="headerlink" href="#ecell4.egfrd.EGFRDFactory.create_simulator" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a EGFRDSimulator instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : EGFRDWorld</p>
<blockquote>
<div><p>A world</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : Model</p>
<blockquote>
<div><p>A simulation model</p>
</div></blockquote>
<p><strong>arg2</strong> : EGFRDWorld</p>
<blockquote>
<div><p>A world</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">EGFRDSimulator:</p>
<blockquote class="last">
<div><p>The created simulator</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDFactory.create_world">
<code class="descname">create_world</code><span class="sig-paren">(</span><em>arg1=None</em><span class="sig-paren">)</span> &rarr; EGFRDWorld<a class="headerlink" href="#ecell4.egfrd.EGFRDFactory.create_world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a EGFRDWorld instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : Real3</p>
<blockquote>
<div><p>The lengths of edges of a EGFRDWorld created</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : str</p>
<blockquote>
<div><p>The path of a HDF5 file for EGFRDWorld</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">EGFRDWorld:</p>
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
<dt id="ecell4.egfrd.EGFRDSimulator">
<em class="property">class </em><code class="descclassname">ecell4.egfrd.</code><code class="descname">EGFRDSimulator</code><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator" title="Permalink to this definition">¶</a></dt>
<dd><p>A class running the simulation with the egfrd algorithm.</p>
<p>EGFRDSimulator(m, w)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.egfrd.EGFRDSimulator.dt">
<code class="descname">dt</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.dt" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the step interval.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDSimulator.initialize">
<code class="descname">initialize</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.initialize" title="Permalink to this definition">¶</a></dt>
<dd><p>Initialize the simulator.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDSimulator.last_reactions">
<code class="descname">last_reactions</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [(ReactionRule, ReactionInfo)]<a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.last_reactions" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDSimulator.model">
<code class="descname">model</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.model" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the model bound.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDSimulator.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the scheduled time for the next step.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDSimulator.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDSimulator.run">
<code class="descname">run</code><span class="sig-paren">(</span><em>duration</em>, <em>observers</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.run" title="Permalink to this definition">¶</a></dt>
<dd><p>Run the simulation.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>duration</strong> : Real</p>
<blockquote>
<div><p>a duration for running a simulation.
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
<dt id="ecell4.egfrd.EGFRDSimulator.set_dt">
<code class="descname">set_dt</code><span class="sig-paren">(</span><em>dt</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.set_dt" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDSimulator.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.set_t" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDSimulator.step">
<code class="descname">step</code><span class="sig-paren">(</span><em>upto=None</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.step" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDSimulator.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDSimulator.world">
<code class="descname">world</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDSimulator.world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the world bound.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.egfrd.EGFRDWorld">
<em class="property">class </em><code class="descclassname">ecell4.egfrd.</code><code class="descname">EGFRDWorld</code><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld" title="Permalink to this definition">¶</a></dt>
<dd><p>A class containing the properties of the egfrd world.</p>
<p>EGFRDWorld(edge_lengths=None, matrix_sizes=None, GSLRandomNumberGenerator rng=None)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.egfrd.EGFRDWorld.add_molecules">
<code class="descname">add_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em>, <em>shape=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.add_molecules" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.apply_boundary">
<code class="descname">apply_boundary</code><span class="sig-paren">(</span><em>Real3 pos</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.apply_boundary" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a position within the world by applying periodic boundaries
to the given position.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDWorld.bind_to">
<code class="descname">bind_to</code><span class="sig-paren">(</span><em>m</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.bind_to" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.distance">
<code class="descname">distance</code><span class="sig-paren">(</span><em>Real3 pos1</em>, <em>Real3 pos2</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.distance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the closest distance between the given positions.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDWorld.edge_lengths">
<code class="descname">edge_lengths</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.edge_lengths" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the edge lengths of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDWorld.get_particle">
<code class="descname">get_particle</code><span class="sig-paren">(</span><em>pid) -&gt; (ParticleID</em>, <em>Particle</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.get_particle" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.has_particle">
<code class="descname">has_particle</code><span class="sig-paren">(</span><em>pid</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.has_particle" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.has_species">
<code class="descname">has_species</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.has_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Check if the given species is in the space or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A species to be found.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>True if the species in the space.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDWorld.list_particles">
<code class="descname">list_particles</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.list_particles" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.list_particles_exact">
<code class="descname">list_particles_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.list_particles_exact" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.list_particles_within_radius">
<code class="descname">list_particles_within_radius</code><span class="sig-paren">(</span><em>pos</em>, <em>radius</em>, <em>ignore1=None</em>, <em>ignore2=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.list_particles_within_radius" title="Permalink to this definition">¶</a></dt>
<dd><p>-&gt; [((ParticleID, Particle), Real)]</p>
<p>Returns a list of pairs of a particle and distance within the given sphere.
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
<dt id="ecell4.egfrd.EGFRDWorld.load">
<code class="descname">load</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.load" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.new_particle">
<code class="descname">new_particle</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2=None) -&gt; (ParticleID</em>, <em>Particle</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.new_particle" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.num_molecules">
<code class="descname">num_molecules</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.num_molecules" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.num_molecules_exact">
<code class="descname">num_molecules_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.num_molecules_exact" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.num_particles">
<code class="descname">num_particles</code><span class="sig-paren">(</span><em>sp=None</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.num_particles" title="Permalink to this definition">¶</a></dt>
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
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first last">Integer: The number of particles (of the given species)</p>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDWorld.num_particles_exact">
<code class="descname">num_particles_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.num_particles_exact" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.remove_molecules">
<code class="descname">remove_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.remove_molecules" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.remove_particle">
<code class="descname">remove_particle</code><span class="sig-paren">(</span><em>pid</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.remove_particle" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.rng">
<code class="descname">rng</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.rng" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a random number generator object.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDWorld.save">
<code class="descname">save</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.save" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.set_t" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time of the world.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.egfrd.EGFRDWorld.update_particle">
<code class="descname">update_particle</code><span class="sig-paren">(</span><em>pid</em>, <em>p</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.update_particle" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.egfrd.EGFRDWorld.volume">
<code class="descname">volume</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.egfrd.EGFRDWorld.volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the volume of the world.</p>
</dd></dl>

</dd></dl>
