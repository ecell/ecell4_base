
<p>Contents:</p>
<div class="toctree-wrapper compound">
<ul class="simple">
</ul>
</div>
<span class="target" id="module-ecell4.core"></span><dl class="class">
<dt id="ecell4.core.ReactionRule">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">ReactionRule</code><a class="headerlink" href="#ecell4.core.ReactionRule" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a reaction rule between <code class="docutils literal"><span class="pre">Species</span></code>.</p>
<p>ReactionRule(reactants=None, products=None, k=None)</p>
<dl class="method">
<dt id="ecell4.core.ReactionRule.add_product">
<code class="descname">add_product</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.add_product" title="Permalink to this definition">¶</a></dt>
<dd><p>Append a product to the end.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>sp (Species): A new product.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.add_reactant">
<code class="descname">add_reactant</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.add_reactant" title="Permalink to this definition">¶</a></dt>
<dd><p>Append a reactant to the end.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>sp (Species): A new reactant.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.as_string">
<code class="descname">as_string</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.as_string" title="Permalink to this definition">¶</a></dt>
<dd><p>Return an unicode string describing this object.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd>str: An unicode string describing this object.</dd>
<dt>Examples:</dt>
<dd><p class="first">The string consists of a list of reactants, a list of products,
and a kinetic rate constant.</p>
<div class="last highlight-python"><div class="highlight"><pre><span class="gp">&gt;&gt;&gt; </span><span class="n">rr</span> <span class="o">=</span> <span class="n">ReactionRule</span><span class="p">([</span><span class="n">Species</span><span class="p">(</span><span class="s">&quot;A&quot;</span><span class="p">),</span> <span class="n">Species</span><span class="p">(</span><span class="s">&quot;B&quot;</span><span class="p">)],</span> <span class="p">[</span><span class="n">Species</span><span class="p">(</span><span class="s">&quot;C&quot;</span><span class="p">)],</span> <span class="mf">1.0</span><span class="p">)</span>
<span class="gp">&gt;&gt;&gt; </span><span class="n">rr</span><span class="o">.</span><span class="n">as_string</span><span class="p">()</span>
<span class="go">u&#39;A+B&gt;C|1&#39;</span>
</pre></div>
</div>
</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.count">
<code class="descname">count</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.count" title="Permalink to this definition">¶</a></dt>
<dd><p>Count the number of matches for reactants.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd><dl class="first last docutils">
<dt>reactants (list): A list of <code class="docutils literal"><span class="pre">Species</span></code>. The order of <code class="docutils literal"><span class="pre">reactants</span></code></dt>
<dd>is respected.</dd>
</dl>
</dd>
<dt>Return:</dt>
<dd>int: The number of matches.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.generate">
<code class="descname">generate</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.generate" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate <a href="#id1"><span class="problematic" id="id2">``</span></a>ReactionRule``s from given reactants.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd><dl class="first last docutils">
<dt>reactants (list): A list of <code class="docutils literal"><span class="pre">Species</span></code>. The order of <code class="docutils literal"><span class="pre">reactants</span></code></dt>
<dd>is respected.</dd>
</dl>
</dd>
<dt>Return:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of <a href="#id3"><span class="problematic" id="id4">``</span></a>ReactionRule``s. The reactants of each</dt>
<dd><code class="docutils literal"><span class="pre">ReactionRule</span></code> are equal to the given <code class="docutils literal"><span class="pre">reactants</span></code>.
If the <code class="docutils literal"><span class="pre">ReactionRule</span></code> does not match the <code class="docutils literal"><span class="pre">reactants</span></code>,
return an empty list.</dd>
</dl>
</dd>
</dl>
<p>Examples:</p>
<div class="highlight-python"><div class="highlight"><pre><span class="gp">&gt;&gt;&gt; </span><span class="n">rr</span> <span class="o">=</span> <span class="n">ReactionRule</span><span class="p">([</span><span class="n">Species</span><span class="p">(</span><span class="s">&quot;_(b=x)&quot;</span><span class="p">)],</span> <span class="p">[</span><span class="n">Species</span><span class="p">(</span><span class="s">&quot;_(b=y)&quot;</span><span class="p">)],</span> <span class="mf">1.0</span><span class="p">)</span>
<span class="gp">&gt;&gt;&gt; </span><span class="n">reactants</span> <span class="o">=</span> <span class="p">[</span><span class="n">Species</span><span class="p">(</span><span class="s">&quot;A(a^1,b=x).B(a^1,b=x)&quot;</span><span class="p">)]</span>
<span class="gp">&gt;&gt;&gt; </span><span class="p">[</span><span class="n">r</span><span class="o">.</span><span class="n">as_string</span><span class="p">()</span> <span class="k">for</span> <span class="n">r</span> <span class="ow">in</span> <span class="n">rr</span><span class="o">.</span><span class="n">generate</span><span class="p">(</span><span class="n">reactants</span><span class="p">)]</span>
<span class="go">[u&#39;A(a^1,b=x).B(a^1,b=x)&gt;A(a^1,b=y).B(a^1,b=x)|1&#39;,</span>
<span class="go"> u&#39;A(a^1,b=x).B(a^1,b=x)&gt;A(a^1,b=x).B(a^1,b=y)|1&#39;]</span>
</pre></div>
</div>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.k">
<code class="descname">k</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.k" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the kinetic rate constant as a float value.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.products">
<code class="descname">products</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.products" title="Permalink to this definition">¶</a></dt>
<dd><p>List all products.</p>
<dl class="docutils">
<dt>Return:</dt>
<dd>list: A list of product <code class="docutils literal"><span class="pre">Species</span></code>.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.reactants">
<code class="descname">reactants</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.reactants" title="Permalink to this definition">¶</a></dt>
<dd><p>List all reactants.</p>
<dl class="docutils">
<dt>Return:</dt>
<dd>list: A list of reactant <code class="docutils literal"><span class="pre">Species</span></code>.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.set_k">
<code class="descname">set_k</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.set_k" title="Permalink to this definition">¶</a></dt>
<dd><p>Set a kinetic rate constant.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>k (float): A kinetic rate constant.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.set_ratelaw">
<code class="descname">set_ratelaw</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.set_ratelaw" title="Permalink to this definition">¶</a></dt>
<dd><p>Warning: This member function will be deprecated.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.set_ratelaw_massaction">
<code class="descname">set_ratelaw_massaction</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.set_ratelaw_massaction" title="Permalink to this definition">¶</a></dt>
<dd><p>Warning: This member function will be deprecated.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Species">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Species</code><a class="headerlink" href="#ecell4.core.Species" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a type of molecules with attributes.</p>
<p>Species(serial=None, radius=None, D=None, location=None)</p>
<dl class="method">
<dt id="ecell4.core.Species.add_unit">
<code class="descname">add_unit</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.add_unit" title="Permalink to this definition">¶</a></dt>
<dd><p>Append an <code class="docutils literal"><span class="pre">UnitSpecies</span></code> to the end.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>usp (UnitSpecies): An <code class="docutils literal"><span class="pre">UnitSpecies</span></code> to be added.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.count">
<code class="descname">count</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.count" title="Permalink to this definition">¶</a></dt>
<dd><p>Count the number of matches for a pattern given as a <code class="docutils literal"><span class="pre">Species</span></code>.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pttrn (Species): A pattern to be count.</dd>
<dt>Returns:</dt>
<dd>int: The number of matches.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.deserialize">
<code class="descname">deserialize</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.deserialize" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the serial. All attributes will be kept.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>serial (string): A new serial as an unicode string.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.get_attribute">
<code class="descname">get_attribute</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.get_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Return an attribute as an unicode string.
If no corresponding attribute is found, raise an error.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>name (str): The name of an attribute.</dd>
<dt>Returns:</dt>
<dd>value (str): The value of the attribute.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.has_attribute">
<code class="descname">has_attribute</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.has_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the attribute exists or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>name (str): The name of an attribute.</dd>
<dt>Returns:</dt>
<dd>bool: True if the attribute exists, False otherwise.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.list_attributes">
<code class="descname">list_attributes</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.list_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>List all attributes.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of pairs of name and value.</dt>
<dd><code class="docutils literal"><span class="pre">name</span></code> and <code class="docutils literal"><span class="pre">value</span></code> are given as unicode strings.</dd>
</dl>
</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.num_units">
<code class="descname">num_units</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.num_units" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of <code class="docutils literal"><span class="pre">UnitSpecies</span></code>.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.remove_attribute">
<code class="descname">remove_attribute</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.remove_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove an attribute.
If no corresponding attribute is found, raise an error.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>name (str): The name of an attribute to be removed.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.serial">
<code class="descname">serial</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.serial" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the serial name as an unicode string.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.set_attribute">
<code class="descname">set_attribute</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.set_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Set an attribute.
If existing already, the attribute will be overwritten.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>name (str): The name of an attribute.
value (str): The value of an attribute.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.units">
<code class="descname">units</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.units" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of all <code class="docutils literal"><span class="pre">UnitSpecies</span></code>.</p>
</dd></dl>

</dd></dl>

<dl class="function">
<dt id="ecell4.core.count_rrmatches">
<code class="descclassname">ecell4.core.</code><code class="descname">count_rrmatches</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.count_rrmatches" title="Permalink to this definition">¶</a></dt>
<dd><p>Count the number of matches for a pattern given as a <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd><p class="first">pttrn (ReactionRule): A pattern.
reactants (list): A list of reactants, <code class="docutils literal"><span class="pre">Species</span></code>.</p>
<blockquote class="last">
<div>The order of reactants is respected.</div></blockquote>
</dd>
<dt>Return:</dt>
<dd>int: The number of matches.</dd>
</dl>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.count_spmatches">
<code class="descclassname">ecell4.core.</code><code class="descname">count_spmatches</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.count_spmatches" title="Permalink to this definition">¶</a></dt>
<dd><p>Count the number of matches for a pattern given as a <code class="docutils literal"><span class="pre">Species</span></code>.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pttrn (Species): A pattern.
sp (Species): A target.</dd>
<dt>Return:</dt>
<dd>int: The number of matches.</dd>
<dt>Note:</dt>
<dd>Use <code class="docutils literal"><span class="pre">Species.count</span></code>.</dd>
</dl>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.create_binding_reaction_rule">
<code class="descclassname">ecell4.core.</code><code class="descname">create_binding_reaction_rule</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.create_binding_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Create a binding <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>reactant1 (Species): One of two reactants.
reactant2 (Species): One of two reactants.
product1 (Species): A product.
k (float): A kinetic parameter.</dd>
<dt>Note:</dt>
<dd>This is equivalent to <code class="docutils literal"><span class="pre">ReactionRule([reactant1,</span> <span class="pre">reactant2],</span> <span class="pre">[product1],</span> <span class="pre">k)</span></code>.</dd>
</dl>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.create_degradation_reaction_rule">
<code class="descclassname">ecell4.core.</code><code class="descname">create_degradation_reaction_rule</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.create_degradation_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Create a degradation <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>reactant1 (Species): A reactant to be degradated.
k (float): A kinetic parameter.</dd>
<dt>Note:</dt>
<dd>This is equivalent to <code class="docutils literal"><span class="pre">ReactionRule([reactant1],</span> <span class="pre">[],</span> <span class="pre">k)</span></code>.</dd>
</dl>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.create_synthesis_reaction_rule">
<code class="descclassname">ecell4.core.</code><code class="descname">create_synthesis_reaction_rule</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.create_synthesis_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Create a synthesis <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>product1 (Species): A product to be synthesized.
k (float): A kinetic parameter.</dd>
<dt>Note:</dt>
<dd>This is equivalent to <code class="docutils literal"><span class="pre">ReactionRule([],</span> <span class="pre">[product1],</span> <span class="pre">k)</span></code>.</dd>
</dl>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.create_unbinding_reaction_rule">
<code class="descclassname">ecell4.core.</code><code class="descname">create_unbinding_reaction_rule</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.create_unbinding_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Create an unbinding <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>reactant1 (Species): A reactant.
product1 (Species): One of two products.
product2 (Species): One of two products.
k (float): A kinetic parameter.</dd>
<dt>Note:</dt>
<dd>This is equivalent to <code class="docutils literal"><span class="pre">ReactionRule([reactant1],</span> <span class="pre">[product1,</span> <span class="pre">product2],</span> <span class="pre">k)</span></code>.</dd>
</dl>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.create_unimolecular_reaction_rule">
<code class="descclassname">ecell4.core.</code><code class="descname">create_unimolecular_reaction_rule</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.create_unimolecular_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Create an unimolecular <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>reactant1 (Species): A reactant to be modified.
product1 (Species): A product.
k (float): A kinetic parameter.</dd>
<dt>Note:</dt>
<dd>This is equivalent to <code class="docutils literal"><span class="pre">ReactionRule([reactant1],</span> <span class="pre">[product1],</span> <span class="pre">k)</span></code>.</dd>
</dl>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.rrgenerate">
<code class="descclassname">ecell4.core.</code><code class="descname">rrgenerate</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.rrgenerate" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate a list of products from the given list of reactants.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd><p class="first">pttrn (ReactionRule): A pattern.
reactants (list): A list of <code class="docutils literal"><span class="pre">Species</span></code>. The order of <code class="docutils literal"><span class="pre">reactants</span></code></p>
<blockquote class="last">
<div>is respected.</div></blockquote>
</dd>
<dt>Return:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of products.</dt>
<dd>The size of the list is equal to the number of matches.
Each element of the list is a list of <code class="docutils literal"><span class="pre">Species</span></code>.</dd>
</dl>
</dd>
<dt>Note:</dt>
<dd>Use <code class="docutils literal"><span class="pre">ReactionRule.generate</span></code>.</dd>
</dl>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.rrmatch">
<code class="descclassname">ecell4.core.</code><code class="descname">rrmatch</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.rrmatch" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if a pattern matches the reactants or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd><p class="first">pttrn (ReactionRule): A pattern.
reactants (list): A list of reactants, <code class="docutils literal"><span class="pre">Species</span></code>.</p>
<blockquote class="last">
<div>The order of reactants is respected.</div></blockquote>
</dd>
<dt>Return:</dt>
<dd><dl class="first last docutils">
<dt>bool: True if <code class="docutils literal"><span class="pre">pttrn</span></code> matches <code class="docutils literal"><span class="pre">reactants</span></code> at least one time,</dt>
<dd>False otherwise.</dd>
</dl>
</dd>
</dl>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.spmatch">
<code class="descclassname">ecell4.core.</code><code class="descname">spmatch</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.spmatch" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if a pattern matches the target <code class="docutils literal"><span class="pre">Species</span></code> or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pttrn (Species): A pattern.
sp (Species): A target.</dd>
<dt>Return:</dt>
<dd>bool: True if <code class="docutils literal"><span class="pre">pttrn</span></code> matches <code class="docutils literal"><span class="pre">sp</span></code> at least one time, False otherwise.</dd>
</dl>
</dd></dl>
