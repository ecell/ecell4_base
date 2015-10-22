  <span class="target" id="module-ecell4.util"></span><dl class="function">
<dt id="ecell4.util.run_simulation">
<code class="descclassname">ecell4.util.</code><code class="descname">run_simulation</code><span class="sig-paren">(</span><em>t</em>, <em>y0={}</em>, <em>volume=1.0</em>, <em>model=None</em>, <em>solver='ode'</em>, <em>factory=None</em>, <em>is_netfree=False</em>, <em>species_list=None</em>, <em>without_reset=False</em>, <em>return_type='matplotlib'</em>, <em>plot_args={}</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.run_simulation" title="Permalink to this definition">¶</a></dt>
<dd><p>Run a simulation with the given model and plot the result on IPython
notebook with matplotlib.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>t</strong> : array</p>
<blockquote>
<div><p>A sequence of time points for which to solve for &#8216;m&#8217;.</p>
</div></blockquote>
<p><strong>y0</strong> : dict</p>
<blockquote>
<div><p>Initial condition.</p>
</div></blockquote>
<p><strong>volume</strong> : Real, optional</p>
<p><strong>model</strong> : Model, optional</p>
<p><strong>solver</strong> : str, optional</p>
<blockquote>
<div><p>Solver type. Choose one from &#8216;ode&#8217;, &#8216;gillespie&#8217;, &#8216;lattice&#8217;, &#8216;meso&#8217;,
&#8216;bd&#8217; and &#8216;egfrd&#8217;. Default is &#8216;ode&#8217;.</p>
</div></blockquote>
<p><strong>species_list</strong> : list of str, optional</p>
<blockquote>
<div><p>A list of names of Species observed. If None, log all.
Default is None.</p>
</div></blockquote>
<p><strong>return_type</strong> : str, optional</p>
<blockquote>
<div><p>Choose a type of return value from &#8216;array&#8217;, &#8216;observer&#8217;,
&#8216;matplotlib&#8217;, &#8216;nyaplot&#8217; or None.
If None, return and plot nothing. Default is &#8216;matplotlib&#8217;.</p>
</div></blockquote>
<p><strong>plot_args</strong> : dict, optional</p>
<blockquote>
<div><p>Arguments for plotting. If plot_type is None, just ignored.</p>
</div></blockquote>
<p><strong>factory</strong> : Factory, optional</p>
<p><strong>is_netfree</strong> : bool, optional</p>
<blockquote>
<div><p>Whether the model is netfree or not. When a model is given as an
argument, just ignored. Default is False.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : list, TimingNumberObserver, or None</p>
<blockquote class="last">
<div><p>Return a value suggested by <code class="docutils literal"><span class="pre">return_type</span></code>.
When <code class="docutils literal"><span class="pre">return_type</span></code> is &#8216;array&#8217;, return a time course data.
When <code class="docutils literal"><span class="pre">return_type</span></code> is &#8216;observer&#8217;, return an observer.
Return nothing if else.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.util.load_world">
<code class="descclassname">ecell4.util.</code><code class="descname">load_world</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.load_world" title="Permalink to this definition">¶</a></dt>
<dd><p>Load a world from the given HDF5 filename.
The return type is determined by <code class="docutils literal"><span class="pre">ecell4.core.load_version_information</span></code>.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>filename</strong> : str</p>
<blockquote>
<div><p>A HDF5 filename.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>w</strong> : World</p>
<blockquote class="last">
<div><p>Return one from <code class="docutils literal"><span class="pre">BDWorld</span></code>, <code class="docutils literal"><span class="pre">EGFRDWorld</span></code>, <code class="docutils literal"><span class="pre">MesoscopicWorld</span></code>,
<code class="docutils literal"><span class="pre">ODEWorld</span></code>, <code class="docutils literal"><span class="pre">GillespieWorld</span></code> and <code class="docutils literal"><span class="pre">LatticeWorld</span></code>.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.util.get_model">
<code class="descclassname">ecell4.util.</code><code class="descname">get_model</code><span class="sig-paren">(</span><em>is_netfree=False</em>, <em>without_reset=False</em>, <em>seeds=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.get_model" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate a model with parameters in the global scope, <code class="docutils literal"><span class="pre">SPECIES_ATTRIBUTES</span></code>
and <code class="docutils literal"><span class="pre">REACTIONRULES</span></code>.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>is_netfree</strong> : bool, optional</p>
<blockquote>
<div><p>Return <code class="docutils literal"><span class="pre">NetfreeModel</span></code> if True, and <code class="docutils literal"><span class="pre">NetworkModel</span></code> if else.
Default is False.</p>
</div></blockquote>
<p><strong>without_reset</strong> : bool, optional</p>
<blockquote>
<div><p>Do not reset the global variables after the generation if True.
Default is False.</p>
</div></blockquote>
<p><strong>seeds</strong> : list, optional</p>
<blockquote>
<div><p>A list of seed <code class="docutils literal"><span class="pre">Species</span></code> for expanding the model.
If this is not None, generate a <code class="docutils literal"><span class="pre">NetfreeModel</span></code> once, and return a
<code class="docutils literal"><span class="pre">NetworkModel</span></code>, which is an expanded form of that with the given seeds.
Default is None.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first last"><strong>model</strong> : NetworkModel, NetfreeModel, or ODENetworkModel</p>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.util.reset_model">
<code class="descclassname">ecell4.util.</code><code class="descname">reset_model</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.reset_model" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset all values, <code class="docutils literal"><span class="pre">SPECIES_ATTRIBUTES</span></code> and <code class="docutils literal"><span class="pre">REACTIONRULES</span></code>,
in the global scope.</p>
</dd></dl>
