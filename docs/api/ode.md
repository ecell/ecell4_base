  <span class="target" id="module-ecell4.ode"></span><dl class="class">
<dt id="ecell4.ode.ODEFactory">
<em class="property">class </em><code class="descclassname">ecell4.ode.</code><code class="descname">ODEFactory</code><a class="headerlink" href="#ecell4.ode.ODEFactory" title="Permalink to this definition">¶</a></dt>
<dd><p>A factory class creating a ODEWorld instance and a ODESimulator instance.</p>
<p>ODEFactory(solvertype=None, dt=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEFactory.create_simulator" title="ecell4.ode.ODEFactory.create_simulator"><code class="xref py py-obj docutils literal"><span class="pre">create_simulator</span></code></a>((arg1,&nbsp;arg2)&nbsp;-&gt;&nbsp;ODESimulator)</td>
<td>Return a ODESimulator instance.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEFactory.create_world" title="ecell4.ode.ODEFactory.create_world"><code class="xref py py-obj docutils literal"><span class="pre">create_world</span></code></a>((arg1=None)&nbsp;-&gt;&nbsp;ODEWorld)</td>
<td>Return a ODEWorld instance.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.ode.ODEFactory.create_simulator">
<code class="descname">create_simulator</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2</em><span class="sig-paren">)</span> &rarr; ODESimulator<a class="headerlink" href="#ecell4.ode.ODEFactory.create_simulator" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a ODESimulator instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : ODEWorld</p>
<blockquote>
<div><p>a world</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : ODENetworkModel or NetworkModel</p>
<blockquote>
<div><p>a simulation model</p>
</div></blockquote>
<p><strong>arg2</strong> : ODEWorld</p>
<blockquote>
<div><p>a world</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">ODESimulator:</p>
<blockquote class="last">
<div><p>the created simulator</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEFactory.create_world">
<code class="descname">create_world</code><span class="sig-paren">(</span><em>arg1=None</em><span class="sig-paren">)</span> &rarr; ODEWorld<a class="headerlink" href="#ecell4.ode.ODEFactory.create_world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a ODEWorld instance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>arg1</strong> : Real3</p>
<blockquote>
<div><p>The lengths of edges of a ODEWorld created</p>
</div></blockquote>
<p><strong>or</strong></p>
<p><strong>arg1</strong> : str</p>
<blockquote>
<div><p>The path of a HDF5 file for ODEWorld</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">ODEWorld:</p>
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
<dt id="ecell4.ode.ODENetworkModel">
<em class="property">class </em><code class="descclassname">ecell4.ode.</code><code class="descname">ODENetworkModel</code><a class="headerlink" href="#ecell4.ode.ODENetworkModel" title="Permalink to this definition">¶</a></dt>
<dd><p>A network model class for ODE simulations.</p>
<p>ODENetworkModel(NetworkModel m=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODENetworkModel.add_reaction_rule" title="ecell4.ode.ODENetworkModel.add_reaction_rule"><code class="xref py py-obj docutils literal"><span class="pre">add_reaction_rule</span></code></a>(rr)</td>
<td>Add a new reaction rule.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODENetworkModel.has_network_model" title="ecell4.ode.ODENetworkModel.has_network_model"><code class="xref py py-obj docutils literal"><span class="pre">has_network_model</span></code></a></td>
<td>Return if this model is bound to a NetworkModel or not.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODENetworkModel.list_species" title="ecell4.ode.ODENetworkModel.list_species"><code class="xref py py-obj docutils literal"><span class="pre">list_species</span></code></a></td>
<td>Return a list of species, contained in reaction rules in the model.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODENetworkModel.num_reaction_rules" title="ecell4.ode.ODENetworkModel.num_reaction_rules"><code class="xref py py-obj docutils literal"><span class="pre">num_reaction_rules</span></code></a></td>
<td>Return a number of reaction rules contained in the model.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODENetworkModel.ode_reaction_rules" title="ecell4.ode.ODENetworkModel.ode_reaction_rules"><code class="xref py py-obj docutils literal"><span class="pre">ode_reaction_rules</span></code></a>(()&nbsp;-&gt;&nbsp;[ODEReactionRule])</td>
<td>Return a list of ODE reaction rules.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODENetworkModel.update_model" title="ecell4.ode.ODENetworkModel.update_model"><code class="xref py py-obj docutils literal"><span class="pre">update_model</span></code></a></td>
<td>Update self to fit the given NetworkModel.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.ode.ODENetworkModel.add_reaction_rule">
<code class="descname">add_reaction_rule</code><span class="sig-paren">(</span><em>rr</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODENetworkModel.add_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a new reaction rule.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>rr</strong> : ReactionRule or ODEReactionRule</p>
<blockquote class="last">
<div><p>A new reaction rule.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODENetworkModel.has_network_model">
<code class="descname">has_network_model</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODENetworkModel.has_network_model" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if this model is bound to a NetworkModel or not.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODENetworkModel.list_species">
<code class="descname">list_species</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODENetworkModel.list_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of species, contained in reaction rules in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODENetworkModel.num_reaction_rules">
<code class="descname">num_reaction_rules</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODENetworkModel.num_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a number of reaction rules contained in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODENetworkModel.ode_reaction_rules">
<code class="descname">ode_reaction_rules</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [ODEReactionRule]<a class="headerlink" href="#ecell4.ode.ODENetworkModel.ode_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of ODE reaction rules.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODENetworkModel.update_model">
<code class="descname">update_model</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODENetworkModel.update_model" title="Permalink to this definition">¶</a></dt>
<dd><p>Update self to fit the given NetworkModel.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.ode.ODERatelaw">
<em class="property">class </em><code class="descclassname">ecell4.ode.</code><code class="descname">ODERatelaw</code><a class="headerlink" href="#ecell4.ode.ODERatelaw" title="Permalink to this definition">¶</a></dt>
<dd><p>An abstract base class for ratelaws bound to ODEReactionRule.</p>
<p>ODERatelaw()</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODERatelaw.as_base" title="ecell4.ode.ODERatelaw.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Return self as a base class.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.ode.ODERatelaw.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODERatelaw.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Return self as a base class. Only for developmental use.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.ode.ODERatelawCallback">
<em class="property">class </em><code class="descclassname">ecell4.ode.</code><code class="descname">ODERatelawCallback</code><a class="headerlink" href="#ecell4.ode.ODERatelawCallback" title="Permalink to this definition">¶</a></dt>
<dd><p>A class for general ratelaws with a callback.</p>
<p>ODERatelawCallback(pyfunc)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODERatelawCallback.as_base" title="ecell4.ode.ODERatelawCallback.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Return self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODERatelawCallback.set_callback" title="ecell4.ode.ODERatelawCallback.set_callback"><code class="xref py py-obj docutils literal"><span class="pre">set_callback</span></code></a>(pyfunc)</td>
<td><table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"></td>
</tr>
</tbody>
</table>
</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.ode.ODERatelawCallback.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODERatelawCallback.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Return self as a base class. Only for developmental use.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODERatelawCallback.set_callback">
<code class="descname">set_callback</code><span class="sig-paren">(</span><em>pyfunc</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODERatelawCallback.set_callback" title="Permalink to this definition">¶</a></dt>
<dd><table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pyfunc</strong> : function</p>
<blockquote class="last">
<div><p>A Python function for the callback
The function must accept five arguments, and return a velocity.
The number of reactants, the number of products, a volume,
the current time, and a ODEReactionRule are given as the
arguments in this order.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Examples</p>
<p>The following callback represents a simple Michaelis-Menten-like
equation:</p>
<div class="highlight-python"><div class="highlight"><pre><span class="gp">&gt;&gt;&gt; </span><span class="n">rl</span> <span class="o">=</span> <span class="n">ODERatelawCallback</span><span class="p">()</span>
<span class="gp">&gt;&gt;&gt; </span><span class="n">rl</span><span class="o">.</span><span class="n">set_callback</span><span class="p">(</span><span class="k">lambda</span> <span class="n">r</span><span class="p">,</span> <span class="n">p</span><span class="p">,</span> <span class="n">v</span><span class="p">,</span> <span class="n">t</span><span class="p">,</span> <span class="n">rr</span><span class="p">:</span> <span class="mf">2.0</span> <span class="o">*</span> <span class="n">r</span><span class="p">[</span><span class="mi">0</span><span class="p">]</span> <span class="o">*</span> <span class="n">r</span><span class="p">[</span><span class="mi">1</span><span class="p">]</span> <span class="o">/</span> <span class="p">(</span><span class="mf">1.0</span> <span class="o">+</span> <span class="n">r</span><span class="p">[</span><span class="mi">1</span><span class="p">]))</span>
</pre></div>
</div>
<p>Here, we expect that the first reactant is an enzyme,
and that the second one is a substrate.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.ode.ODERatelawMassAction">
<em class="property">class </em><code class="descclassname">ecell4.ode.</code><code class="descname">ODERatelawMassAction</code><a class="headerlink" href="#ecell4.ode.ODERatelawMassAction" title="Permalink to this definition">¶</a></dt>
<dd><p>A class for mass action ratelaws.</p>
<p>ODERatelawMassAction(Real k)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODERatelawMassAction.as_base" title="ecell4.ode.ODERatelawMassAction.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Return self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODERatelawMassAction.get_k" title="ecell4.ode.ODERatelawMassAction.get_k"><code class="xref py py-obj docutils literal"><span class="pre">get_k</span></code></a></td>
<td>Return the kinetic rate constant as a float value.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODERatelawMassAction.is_available" title="ecell4.ode.ODERatelawMassAction.is_available"><code class="xref py py-obj docutils literal"><span class="pre">is_available</span></code></a></td>
<td>Check if this ratelaw is available or not.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODERatelawMassAction.set_k" title="ecell4.ode.ODERatelawMassAction.set_k"><code class="xref py py-obj docutils literal"><span class="pre">set_k</span></code></a>(k)</td>
<td>Set a kinetic rate constant.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.ode.ODERatelawMassAction.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODERatelawMassAction.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Return self as a base class. Only for developmental use.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODERatelawMassAction.get_k">
<code class="descname">get_k</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODERatelawMassAction.get_k" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the kinetic rate constant as a float value.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODERatelawMassAction.is_available">
<code class="descname">is_available</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODERatelawMassAction.is_available" title="Permalink to this definition">¶</a></dt>
<dd><p>Check if this ratelaw is available or not. Return True always.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODERatelawMassAction.set_k">
<code class="descname">set_k</code><span class="sig-paren">(</span><em>k</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODERatelawMassAction.set_k" title="Permalink to this definition">¶</a></dt>
<dd><p>Set a kinetic rate constant.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>k</strong> : float</p>
<blockquote class="last">
<div><p>A kinetic rate constant.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.ode.ODEReactionRule">
<em class="property">class </em><code class="descclassname">ecell4.ode.</code><code class="descname">ODEReactionRule</code><a class="headerlink" href="#ecell4.ode.ODEReactionRule" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a reaction rule between <code class="docutils literal"><span class="pre">Species</span></code>, which accepts at most
one rate law to calculate the flux.</p>
<p>ODEReactionRule()</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.add_product" title="ecell4.ode.ODEReactionRule.add_product"><code class="xref py py-obj docutils literal"><span class="pre">add_product</span></code></a>(sp[,&nbsp;coeff])</td>
<td>Append a product to the end.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.add_reactant" title="ecell4.ode.ODEReactionRule.add_reactant"><code class="xref py py-obj docutils literal"><span class="pre">add_reactant</span></code></a>(sp[,&nbsp;coeff])</td>
<td>Append a reactant to the end.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.as_string" title="ecell4.ode.ODEReactionRule.as_string"><code class="xref py py-obj docutils literal"><span class="pre">as_string</span></code></a>(()&nbsp;-&gt;&nbsp;str)</td>
<td>Return an unicode string describing this object.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.has_ratelaw" title="ecell4.ode.ODEReactionRule.has_ratelaw"><code class="xref py py-obj docutils literal"><span class="pre">has_ratelaw</span></code></a></td>
<td>Return if a ratelaw is bound or not.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.is_massaction" title="ecell4.ode.ODEReactionRule.is_massaction"><code class="xref py py-obj docutils literal"><span class="pre">is_massaction</span></code></a></td>
<td>Return if a mass action ratelaw is bound or not.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.k" title="ecell4.ode.ODEReactionRule.k"><code class="xref py py-obj docutils literal"><span class="pre">k</span></code></a></td>
<td>Return the kinetic rate constant as a float value.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.products" title="ecell4.ode.ODEReactionRule.products"><code class="xref py py-obj docutils literal"><span class="pre">products</span></code></a></td>
<td>List all products.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.products_coefficients" title="ecell4.ode.ODEReactionRule.products_coefficients"><code class="xref py py-obj docutils literal"><span class="pre">products_coefficients</span></code></a>(()&nbsp;-&gt;&nbsp;[Integer])</td>
<td>List all coefficients for products.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.reactants" title="ecell4.ode.ODEReactionRule.reactants"><code class="xref py py-obj docutils literal"><span class="pre">reactants</span></code></a></td>
<td>List all reactants.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.reactants_coefficients" title="ecell4.ode.ODEReactionRule.reactants_coefficients"><code class="xref py py-obj docutils literal"><span class="pre">reactants_coefficients</span></code></a>(()&nbsp;-&gt;&nbsp;[Integer])</td>
<td>List all coefficients for reactants.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.set_k" title="ecell4.ode.ODEReactionRule.set_k"><code class="xref py py-obj docutils literal"><span class="pre">set_k</span></code></a>(k)</td>
<td>Set a kinetic rate constant.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.set_product_coefficient" title="ecell4.ode.ODEReactionRule.set_product_coefficient"><code class="xref py py-obj docutils literal"><span class="pre">set_product_coefficient</span></code></a>(index,&nbsp;coeff)</td>
<td>Set a stoichiometry coefficient of a product at the given index.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.set_ratelaw" title="ecell4.ode.ODEReactionRule.set_ratelaw"><code class="xref py py-obj docutils literal"><span class="pre">set_ratelaw</span></code></a>(ratelaw_obj)</td>
<td>Bind a ratelaw.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.set_ratelaw_massaction" title="ecell4.ode.ODEReactionRule.set_ratelaw_massaction"><code class="xref py py-obj docutils literal"><span class="pre">set_ratelaw_massaction</span></code></a>(ratelaw_obj)</td>
<td>Bind a mass action ratelaw.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEReactionRule.set_reactant_coefficient" title="ecell4.ode.ODEReactionRule.set_reactant_coefficient"><code class="xref py py-obj docutils literal"><span class="pre">set_reactant_coefficient</span></code></a>(index,&nbsp;coeff)</td>
<td>Set a stoichiometry coefficient of a reactant at the given index.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.add_product">
<code class="descname">add_product</code><span class="sig-paren">(</span><em>sp</em>, <em>coeff=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.add_product" title="Permalink to this definition">¶</a></dt>
<dd><p>Append a product to the end.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A new product.</p>
</div></blockquote>
<p><strong>coeff</strong> : Integer</p>
<blockquote class="last">
<div><p>A stoichiometry coefficient.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.add_reactant">
<code class="descname">add_reactant</code><span class="sig-paren">(</span><em>sp</em>, <em>coeff=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.add_reactant" title="Permalink to this definition">¶</a></dt>
<dd><p>Append a reactant to the end.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A new reactant.</p>
</div></blockquote>
<p><strong>coeff</strong> : Integer</p>
<blockquote class="last">
<div><p>A stoichiometry coefficient.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.as_string">
<code class="descname">as_string</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; str<a class="headerlink" href="#ecell4.ode.ODEReactionRule.as_string" title="Permalink to this definition">¶</a></dt>
<dd><p>Return an unicode string describing this object.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">str:</p>
<blockquote class="last">
<div><p>An unicode string describing this object.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.has_ratelaw">
<code class="descname">has_ratelaw</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.has_ratelaw" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if a ratelaw is bound or not.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.is_massaction">
<code class="descname">is_massaction</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.is_massaction" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if a mass action ratelaw is bound or not.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.k">
<code class="descname">k</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.k" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the kinetic rate constant as a float value.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.products">
<code class="descname">products</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.products" title="Permalink to this definition">¶</a></dt>
<dd><p>List all products.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of product <code class="docutils literal"><span class="pre">Species</span></code>.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.products_coefficients">
<code class="descname">products_coefficients</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [Integer]<a class="headerlink" href="#ecell4.ode.ODEReactionRule.products_coefficients" title="Permalink to this definition">¶</a></dt>
<dd><p>List all coefficients for products.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of product coefficients.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.reactants">
<code class="descname">reactants</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.reactants" title="Permalink to this definition">¶</a></dt>
<dd><p>List all reactants.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of reactant <code class="docutils literal"><span class="pre">Species</span></code>.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.reactants_coefficients">
<code class="descname">reactants_coefficients</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [Integer]<a class="headerlink" href="#ecell4.ode.ODEReactionRule.reactants_coefficients" title="Permalink to this definition">¶</a></dt>
<dd><p>List all coefficients for reactants.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of reactant coefficients.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.set_k">
<code class="descname">set_k</code><span class="sig-paren">(</span><em>k</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.set_k" title="Permalink to this definition">¶</a></dt>
<dd><p>Set a kinetic rate constant.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>k</strong> : float</p>
<blockquote class="last">
<div><p>A kinetic rate constant.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.set_product_coefficient">
<code class="descname">set_product_coefficient</code><span class="sig-paren">(</span><em>index</em>, <em>coeff</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.set_product_coefficient" title="Permalink to this definition">¶</a></dt>
<dd><p>Set a stoichiometry coefficient of a product at the given index.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>index</strong> : Integer</p>
<blockquote>
<div><p>An index pointing the target product.</p>
</div></blockquote>
<p><strong>coeff</strong> : Integer</p>
<blockquote class="last">
<div><p>A stoichiometry coefficient.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.set_ratelaw">
<code class="descname">set_ratelaw</code><span class="sig-paren">(</span><em>ratelaw_obj</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.set_ratelaw" title="Permalink to this definition">¶</a></dt>
<dd><p>Bind a ratelaw.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>ratelaw_obj</strong> : ODERatelaw</p>
<blockquote class="last">
<div><p>A ratelaw</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.set_ratelaw_massaction">
<code class="descname">set_ratelaw_massaction</code><span class="sig-paren">(</span><em>ratelaw_obj</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.set_ratelaw_massaction" title="Permalink to this definition">¶</a></dt>
<dd><p>Bind a mass action ratelaw. This will be deprecated soon.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>ratelaw_obj</strong> : ODERatelawMassAction</p>
<blockquote class="last">
<div><p>A ratelaw</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEReactionRule.set_reactant_coefficient">
<code class="descname">set_reactant_coefficient</code><span class="sig-paren">(</span><em>index</em>, <em>coeff</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEReactionRule.set_reactant_coefficient" title="Permalink to this definition">¶</a></dt>
<dd><p>Set a stoichiometry coefficient of a reactant at the given index.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>index</strong> : Integer</p>
<blockquote>
<div><p>An index pointing the target reactant.</p>
</div></blockquote>
<p><strong>coeff</strong> : Integer</p>
<blockquote class="last">
<div><p>A stoichiometry coefficient.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.ode.ODESimulator">
<em class="property">class </em><code class="descclassname">ecell4.ode.</code><code class="descname">ODESimulator</code><a class="headerlink" href="#ecell4.ode.ODESimulator" title="Permalink to this definition">¶</a></dt>
<dd><p>A class running the simulation with the ode algorithm.</p>
<p>ODESimulator(m, w, solver_type)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.absolute_tolerance" title="ecell4.ode.ODESimulator.absolute_tolerance"><code class="xref py py-obj docutils literal"><span class="pre">absolute_tolerance</span></code></a></td>
<td>Return the absolute tolerance.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.dt" title="ecell4.ode.ODESimulator.dt"><code class="xref py py-obj docutils literal"><span class="pre">dt</span></code></a></td>
<td>Return the step interval.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.initialize" title="ecell4.ode.ODESimulator.initialize"><code class="xref py py-obj docutils literal"><span class="pre">initialize</span></code></a></td>
<td>Initialize the simulator.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.next_time" title="ecell4.ode.ODESimulator.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the scheduled time for the next step.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.num_steps" title="ecell4.ode.ODESimulator.num_steps"><code class="xref py py-obj docutils literal"><span class="pre">num_steps</span></code></a></td>
<td>Return the number of steps.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.relative_tolerance" title="ecell4.ode.ODESimulator.relative_tolerance"><code class="xref py py-obj docutils literal"><span class="pre">relative_tolerance</span></code></a></td>
<td>Return the relative tolerance.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.run" title="ecell4.ode.ODESimulator.run"><code class="xref py py-obj docutils literal"><span class="pre">run</span></code></a>(duration,&nbsp;observers)</td>
<td>Run the simulation.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.set_absolute_tolerance" title="ecell4.ode.ODESimulator.set_absolute_tolerance"><code class="xref py py-obj docutils literal"><span class="pre">set_absolute_tolerance</span></code></a>(abs_tol)</td>
<td>Set the absolute tolerance.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.set_dt" title="ecell4.ode.ODESimulator.set_dt"><code class="xref py py-obj docutils literal"><span class="pre">set_dt</span></code></a>(dt)</td>
<td>Set a step interval.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.set_relative_tolerance" title="ecell4.ode.ODESimulator.set_relative_tolerance"><code class="xref py py-obj docutils literal"><span class="pre">set_relative_tolerance</span></code></a>(rel_tol)</td>
<td>Set the relative tolerance.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.set_t" title="ecell4.ode.ODESimulator.set_t"><code class="xref py py-obj docutils literal"><span class="pre">set_t</span></code></a>(t)</td>
<td>Set the current time.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.step" title="ecell4.ode.ODESimulator.step"><code class="xref py py-obj docutils literal"><span class="pre">step</span></code></a>((upto=None)&nbsp;-&gt;&nbsp;bool)</td>
<td>Step the simulation.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODESimulator.t" title="ecell4.ode.ODESimulator.t"><code class="xref py py-obj docutils literal"><span class="pre">t</span></code></a></td>
<td>Return the time.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.ode.ODESimulator.absolute_tolerance">
<code class="descname">absolute_tolerance</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.absolute_tolerance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the absolute tolerance.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODESimulator.dt">
<code class="descname">dt</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.dt" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the step interval.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODESimulator.initialize">
<code class="descname">initialize</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.initialize" title="Permalink to this definition">¶</a></dt>
<dd><p>Initialize the simulator.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODESimulator.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the scheduled time for the next step.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODESimulator.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODESimulator.relative_tolerance">
<code class="descname">relative_tolerance</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.relative_tolerance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the relative tolerance.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODESimulator.run">
<code class="descname">run</code><span class="sig-paren">(</span><em>duration</em>, <em>observers</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.run" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.ode.ODESimulator.set_absolute_tolerance">
<code class="descname">set_absolute_tolerance</code><span class="sig-paren">(</span><em>abs_tol</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.set_absolute_tolerance" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the absolute tolerance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>abs_tol</strong> : Real</p>
<blockquote class="last">
<div><p>an absolute tolerance.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODESimulator.set_dt">
<code class="descname">set_dt</code><span class="sig-paren">(</span><em>dt</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.set_dt" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.ode.ODESimulator.set_relative_tolerance">
<code class="descname">set_relative_tolerance</code><span class="sig-paren">(</span><em>rel_tol</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.set_relative_tolerance" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the relative tolerance.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>rel_tol</strong> : Real</p>
<blockquote class="last">
<div><p>an relative tolerance.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODESimulator.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.set_t" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.ode.ODESimulator.step">
<code class="descname">step</code><span class="sig-paren">(</span><em>upto=None</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.ode.ODESimulator.step" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.ode.ODESimulator.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODESimulator.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the time.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.ode.ODEWorld">
<em class="property">class </em><code class="descclassname">ecell4.ode.</code><code class="descname">ODEWorld</code><a class="headerlink" href="#ecell4.ode.ODEWorld" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing the World for ODE simulations.</p>
<p>ODEWorld(edge_lengths=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.add_molecules" title="ecell4.ode.ODEWorld.add_molecules"><code class="xref py py-obj docutils literal"><span class="pre">add_molecules</span></code></a>(sp,&nbsp;num[,&nbsp;shape])</td>
<td>Add some molecules.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.as_base" title="ecell4.ode.ODEWorld.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Return self as a base class.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.bind_to" title="ecell4.ode.ODEWorld.bind_to"><code class="xref py py-obj docutils literal"><span class="pre">bind_to</span></code></a>(m)</td>
<td>Bind a model.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.edge_lengths" title="ecell4.ode.ODEWorld.edge_lengths"><code class="xref py py-obj docutils literal"><span class="pre">edge_lengths</span></code></a>(()&nbsp;-&gt;&nbsp;Real3)</td>
<td>Return edge lengths for the space.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.get_value" title="ecell4.ode.ODEWorld.get_value"><code class="xref py py-obj docutils literal"><span class="pre">get_value</span></code></a>((sp)&nbsp;-&gt;&nbsp;Real)</td>
<td>Return the value matched to a given species.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.get_value_exact" title="ecell4.ode.ODEWorld.get_value_exact"><code class="xref py py-obj docutils literal"><span class="pre">get_value_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;Real)</td>
<td>Return the value connected to a given species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.has_species" title="ecell4.ode.ODEWorld.has_species"><code class="xref py py-obj docutils literal"><span class="pre">has_species</span></code></a>((sp)&nbsp;-&gt;&nbsp;bool)</td>
<td>Check if the given species is belonging to this.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.list_species" title="ecell4.ode.ODEWorld.list_species"><code class="xref py py-obj docutils literal"><span class="pre">list_species</span></code></a></td>
<td>Return a list of species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.load" title="ecell4.ode.ODEWorld.load"><code class="xref py py-obj docutils literal"><span class="pre">load</span></code></a>(filename)</td>
<td>Load a HDF5 file to the current state.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.num_molecules" title="ecell4.ode.ODEWorld.num_molecules"><code class="xref py py-obj docutils literal"><span class="pre">num_molecules</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of molecules.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.num_molecules_exact" title="ecell4.ode.ODEWorld.num_molecules_exact"><code class="xref py py-obj docutils literal"><span class="pre">num_molecules_exact</span></code></a>((sp)&nbsp;-&gt;&nbsp;Integer)</td>
<td>Return the number of molecules of a given species.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.release_species" title="ecell4.ode.ODEWorld.release_species"><code class="xref py py-obj docutils literal"><span class="pre">release_species</span></code></a>(sp)</td>
<td>Release a value for the given species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.remove_molecules" title="ecell4.ode.ODEWorld.remove_molecules"><code class="xref py py-obj docutils literal"><span class="pre">remove_molecules</span></code></a>(sp,&nbsp;num)</td>
<td>Remove molecules</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.reserve_species" title="ecell4.ode.ODEWorld.reserve_species"><code class="xref py py-obj docutils literal"><span class="pre">reserve_species</span></code></a>(sp)</td>
<td>Reserve a value for the given species.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.save" title="ecell4.ode.ODEWorld.save"><code class="xref py py-obj docutils literal"><span class="pre">save</span></code></a>(filename)</td>
<td>Save the current state to a HDF5 file.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.set_t" title="ecell4.ode.ODEWorld.set_t"><code class="xref py py-obj docutils literal"><span class="pre">set_t</span></code></a>(t)</td>
<td>Set the current time.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.set_value" title="ecell4.ode.ODEWorld.set_value"><code class="xref py py-obj docutils literal"><span class="pre">set_value</span></code></a>(sp,&nbsp;value)</td>
<td>Set the value of the given species.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.set_volume" title="ecell4.ode.ODEWorld.set_volume"><code class="xref py py-obj docutils literal"><span class="pre">set_volume</span></code></a>(volume)</td>
<td>Set a volume.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.t" title="ecell4.ode.ODEWorld.t"><code class="xref py py-obj docutils literal"><span class="pre">t</span></code></a></td>
<td>Return the current time.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.ode.ODEWorld.volume" title="ecell4.ode.ODEWorld.volume"><code class="xref py py-obj docutils literal"><span class="pre">volume</span></code></a></td>
<td>Return a volume.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.ode.ODEWorld.add_molecules">
<code class="descname">add_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em>, <em>shape=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.add_molecules" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.ode.ODEWorld.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Return self as a base class. Only for developmental use.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.bind_to">
<code class="descname">bind_to</code><span class="sig-paren">(</span><em>m</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.bind_to" title="Permalink to this definition">¶</a></dt>
<dd><p>Bind a model.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>m</strong> : ODENetworkModel or NetworkModel</p>
<blockquote class="last">
<div><p>a model to be bound</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.edge_lengths">
<code class="descname">edge_lengths</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.ode.ODEWorld.edge_lengths" title="Permalink to this definition">¶</a></dt>
<dd><p>Return edge lengths for the space.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.get_value">
<code class="descname">get_value</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.ode.ODEWorld.get_value" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the value matched to a given species.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>a pattern whose value you get</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real:</p>
<blockquote class="last">
<div><p>the value matched to a given species</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.get_value_exact">
<code class="descname">get_value_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.ode.ODEWorld.get_value_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the value connected to a given species.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>a species whose value you get</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real:</p>
<blockquote class="last">
<div><p>the value connected to a given species</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.has_species">
<code class="descname">has_species</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.ode.ODEWorld.has_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Check if the given species is belonging to this.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>a species to be checked.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>True if the given species is contained.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.list_species">
<code class="descname">list_species</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.list_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of species.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.load">
<code class="descname">load</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.load" title="Permalink to this definition">¶</a></dt>
<dd><p>Load a HDF5 file to the current state.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>filename</strong> : str</p>
<blockquote class="last">
<div><p>a file name to be loaded.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.num_molecules">
<code class="descname">num_molecules</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.ode.ODEWorld.num_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of molecules. A value is rounded to an integer.
See set_value also.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species, optional</p>
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
<dt id="ecell4.ode.ODEWorld.num_molecules_exact">
<code class="descname">num_molecules_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.ode.ODEWorld.num_molecules_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of molecules of a given species.
A value is rounded to an integer. See get_value_exact also.</p>
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
<dt id="ecell4.ode.ODEWorld.release_species">
<code class="descname">release_species</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.release_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Release a value for the given species.
This function is mainly for developers.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote class="last">
<div><p>a species to be released.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.remove_molecules">
<code class="descname">remove_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.remove_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove molecules</p>
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
<dt id="ecell4.ode.ODEWorld.reserve_species">
<code class="descname">reserve_species</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.reserve_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Reserve a value for the given species. Use set_value.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote class="last">
<div><p>a species to be reserved.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.save">
<code class="descname">save</code><span class="sig-paren">(</span><em>filename</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.save" title="Permalink to this definition">¶</a></dt>
<dd><p>Save the current state to a HDF5 file.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>filename</strong> : str</p>
<blockquote class="last">
<div><p>a file name to be saved.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.set_t" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the current time.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.set_value">
<code class="descname">set_value</code><span class="sig-paren">(</span><em>sp</em>, <em>value</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.set_value" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the value of the given species.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>a species whose value you set</p>
</div></blockquote>
<p><strong>value</strong> : Real</p>
<blockquote class="last">
<div><p>a value set</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.set_volume">
<code class="descname">set_volume</code><span class="sig-paren">(</span><em>volume</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.set_volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Set a volume.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.t">
<code class="descname">t</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.t" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the current time.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.ode.ODEWorld.volume">
<code class="descname">volume</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.ode.ODEWorld.volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a volume.</p>
</dd></dl>

</dd></dl>
