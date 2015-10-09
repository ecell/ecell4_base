  <span class="target" id="module-ecell4.util.decorator"></span><dl class="function">
<dt id="ecell4.util.decorator.get_model">
<code class="descclassname">ecell4.util.decorator.</code><code class="descname">get_model</code><span class="sig-paren">(</span><em>is_netfree=False</em>, <em>without_reset=False</em>, <em>seeds=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.decorator.get_model" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.util.decorator.reset_model">
<code class="descclassname">ecell4.util.decorator.</code><code class="descname">reset_model</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.util.decorator.reset_model" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset all values, <code class="docutils literal"><span class="pre">SPECIES_ATTRIBUTES</span></code> and <code class="docutils literal"><span class="pre">REACTIONRULES</span></code>,
in the global scope.</p>
</dd></dl>
