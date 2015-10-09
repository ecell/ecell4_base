  <span class="target" id="module-ecell4.core"></span><dl class="class">
<dt id="ecell4.core.AABB">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">AABB</code><a class="headerlink" href="#ecell4.core.AABB" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing an axis aligned bounding box (AABB),
which is available to define structures.</p>
<p>AABB(lower, upper)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.AABB.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.AABB.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.AABB.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.AABB.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.AABB.distance">
<code class="descname">distance</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.AABB.distance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a minimum distance from the given point to the surface.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>distance</strong> : float</p>
<blockquote class="last">
<div><p>A minimum distance from the given point.
Negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.AABB.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.AABB.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : float</p>
<blockquote class="last">
<div><p>Zero or negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.AABB.lower">
<code class="descname">lower</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.AABB.lower" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a vertex suggesting the lower bounds.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.AABB.upper">
<code class="descname">upper</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.AABB.upper" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a vertex suggesting the upper bounds.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Cylinder">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Cylinder</code><a class="headerlink" href="#ecell4.core.Cylinder" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a cylinder shape, which is available to define
structures.</p>
<p>Cylinder(center, radius, axis, half_height)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.Cylinder.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Cylinder.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Cylinder.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Cylinder.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Cylinder.distance">
<code class="descname">distance</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Cylinder.distance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a minimum distance from the given point to the surface.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>distance</strong> : float</p>
<blockquote class="last">
<div><p>A minimum distance from the given point.
Negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Cylinder.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Cylinder.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : float</p>
<blockquote class="last">
<div><p>Zero or negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Cylinder.surface">
<code class="descname">surface</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Cylinder.surface" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a surface shape.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>shape</strong> : CylindricalSurface</p>
<blockquote class="last">
<div><p>The surface shape.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.CylindricalSurface">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">CylindricalSurface</code><a class="headerlink" href="#ecell4.core.CylindricalSurface" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a hollow cylindrical surface, which is
available to define structures.</p>
<p>CylindricalSurface(center, radius, axis, half_height)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.CylindricalSurface.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.CylindricalSurface.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.CylindricalSurface.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.CylindricalSurface.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.CylindricalSurface.distance">
<code class="descname">distance</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.CylindricalSurface.distance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a minimum distance from the given point to the surface.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>distance</strong> : float</p>
<blockquote class="last">
<div><p>A minimum distance from the given point.
Negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.CylindricalSurface.inside">
<code class="descname">inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.CylindricalSurface.inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a volume shape.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>shape</strong> : Cylinder</p>
<blockquote class="last">
<div><p>The volume shape.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.CylindricalSurface.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.CylindricalSurface.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : float</p>
<blockquote class="last">
<div><p>Zero or negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.FixedIntervalCSVObserver">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">FixedIntervalCSVObserver</code><a class="headerlink" href="#ecell4.core.FixedIntervalCSVObserver" title="Permalink to this definition">¶</a></dt>
<dd><p>An <code class="docutils literal"><span class="pre">Observer</span></code> class to log the state of <code class="docutils literal"><span class="pre">World</span></code> in CSV format
with the fixed step interval.
This <code class="docutils literal"><span class="pre">Observer</span></code> saves the <code class="docutils literal"><span class="pre">World</span></code> at the current time first, and
then keeps saving every after the interval.</p>
<p>FixedIntervalCSVObserver(dt, filename)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.FixedIntervalCSVObserver.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalCSVObserver.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalCSVObserver.filename">
<code class="descname">filename</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalCSVObserver.filename" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a file name to be saved at the next time</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalCSVObserver.log">
<code class="descname">log</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalCSVObserver.log" title="Permalink to this definition">¶</a></dt>
<dd><p>Force to log the given <code class="docutils literal"><span class="pre">World</span></code> to a file.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>w</strong> : Space</p>
<blockquote class="last">
<div><p>A <code class="docutils literal"><span class="pre">Space</span></code> (<code class="docutils literal"><span class="pre">World</span></code>) to be logged.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Examples</p>
<p>This is an easy way to save a <code class="docutils literal"><span class="pre">World</span></code> in CSV format without
running a simulation.</p>
<div class="highlight-python"><div class="highlight"><pre><span class="gp">&gt;&gt;&gt; </span><span class="n">w</span> <span class="o">=</span> <span class="n">lattice</span><span class="o">.</span><span class="n">LatticeWorld</span><span class="p">(</span><span class="n">Real3</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">),</span> <span class="mf">0.005</span><span class="p">)</span>
<span class="gp">&gt;&gt;&gt; </span><span class="n">w</span><span class="o">.</span><span class="n">bind_to</span><span class="p">(</span><span class="n">NetworkModel</span><span class="p">())</span>
<span class="gp">&gt;&gt;&gt; </span><span class="n">w</span><span class="o">.</span><span class="n">add_molecules</span><span class="p">(</span><span class="n">Species</span><span class="p">(</span><span class="s">&quot;A&quot;</span><span class="p">),</span> <span class="mi">3</span><span class="p">)</span>
<span class="gp">&gt;&gt;&gt; </span><span class="n">FixedIntervalCSVObserver</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="s">&quot;test.csv&quot;</span><span class="p">)</span><span class="o">.</span><span class="n">log</span><span class="p">(</span><span class="n">w</span><span class="p">)</span>
<span class="gp">&gt;&gt;&gt; </span><span class="k">print</span><span class="p">(</span><span class="nb">open</span><span class="p">(</span><span class="s">&quot;test.csv&quot;</span><span class="p">)</span><span class="o">.</span><span class="n">read</span><span class="p">())</span>
<span class="go">x,y,z,r,sid</span>
<span class="go">0.10614455552060439,0.66106605822212161,0.81500000000000006,0.0050000000000000001,0</span>
<span class="go">0.38375339303603129,0.37527767497325676,0.23999999999999999,0.0050000000000000001,0</span>
<span class="go">0.25311394008759508,0.05484827557301445,0.495,0.0050000000000000001,0</span>
</pre></div>
</div>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalCSVObserver.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalCSVObserver.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the next time for logging.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalCSVObserver.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalCSVObserver.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalCSVObserver.reset">
<code class="descname">reset</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalCSVObserver.reset" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the internal state.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.FixedIntervalHDF5Observer">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">FixedIntervalHDF5Observer</code><a class="headerlink" href="#ecell4.core.FixedIntervalHDF5Observer" title="Permalink to this definition">¶</a></dt>
<dd><p>An <code class="docutils literal"><span class="pre">Observer</span></code> class to log the state of <code class="docutils literal"><span class="pre">World</span></code> in HDF5 format
with the fixed step interval.
This <code class="docutils literal"><span class="pre">Observer</span></code> saves the <code class="docutils literal"><span class="pre">World</span></code> at the current time first, and
then keeps saving every after the interval.</p>
<p>FixedIntervalHDF5Observer(dt, filename)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.FixedIntervalHDF5Observer.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalHDF5Observer.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalHDF5Observer.filename">
<code class="descname">filename</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalHDF5Observer.filename" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a file name to be saved at the next time</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalHDF5Observer.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalHDF5Observer.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the next time for logging.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalHDF5Observer.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalHDF5Observer.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalHDF5Observer.reset">
<code class="descname">reset</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalHDF5Observer.reset" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the internal state.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.FixedIntervalNumberObserver">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">FixedIntervalNumberObserver</code><a class="headerlink" href="#ecell4.core.FixedIntervalNumberObserver" title="Permalink to this definition">¶</a></dt>
<dd><p>An <code class="docutils literal"><span class="pre">Observer``class</span> <span class="pre">to</span> <span class="pre">log</span> <span class="pre">the</span> <span class="pre">number</span> <span class="pre">of</span> <span class="pre">molecules</span> <span class="pre">with</span> <span class="pre">the</span> <span class="pre">fixed</span>
<span class="pre">step</span> <span class="pre">interval.</span>
<span class="pre">This</span> <span class="pre">``Observer</span></code> logs at the current time first, and then keeps logging
every after the interval.</p>
<p>FixedIntervalNumberObserver(dt, species)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.FixedIntervalNumberObserver.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalNumberObserver.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalNumberObserver.data">
<code class="descname">data</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalNumberObserver.data" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of the number of molecules you specified.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of lists of the numbers of molecules.
The size of a return value is equal to <code class="docutils literal"><span class="pre">num_steps</span></code>.
Each element of a return value is a list consisting of
time and the number of molecules specified at the construction.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalNumberObserver.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalNumberObserver.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the next time for logging.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalNumberObserver.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalNumberObserver.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalNumberObserver.reset">
<code class="descname">reset</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalNumberObserver.reset" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the internal state.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalNumberObserver.targets">
<code class="descname">targets</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalNumberObserver.targets" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of <code class="docutils literal"><span class="pre">Species</span></code>, which this <code class="docutils literal"><span class="pre">Observer</span></code> observes</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of <code class="docutils literal"><span class="pre">Species</span></code>. This is generated from arguments
you gave at the construction.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.FixedIntervalTrajectoryObserver">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">FixedIntervalTrajectoryObserver</code><a class="headerlink" href="#ecell4.core.FixedIntervalTrajectoryObserver" title="Permalink to this definition">¶</a></dt>
<dd><p>An <code class="docutils literal"><span class="pre">Observer</span></code> class to trace and log trajectories of diffusing
particles in a <code class="docutils literal"><span class="pre">World</span></code> with the fixed step interval.
This <code class="docutils literal"><span class="pre">Observer</span></code> logs at the current time first, and then keeps logging
every after the interval.</p>
<p>FixedIntervalTrajectoryObserver(dt, pids, resolve_boundary=None)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.FixedIntervalTrajectoryObserver.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalTrajectoryObserver.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalTrajectoryObserver.data">
<code class="descname">data</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalTrajectoryObserver.data" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of trajectories for each particles.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of lists of <code class="docutils literal"><span class="pre">Real3</span></code>. An element of a return value
is corresponding the trajectory of each particle. Thus, the size
of a return value is the same with that of <code class="docutils literal"><span class="pre">pids</span></code> you gave
at the construction.
If a particle corresponding to the given <code class="docutils literal"><span class="pre">ParticleID</span></code> is missing,
i.e. for a reaction, this <code class="docutils literal"><span class="pre">Observer</span></code> just skips to log the
position. Therefore, lengths of the trajectories can be diverse.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalTrajectoryObserver.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalTrajectoryObserver.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the next time for logging.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalTrajectoryObserver.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalTrajectoryObserver.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalTrajectoryObserver.reset">
<code class="descname">reset</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalTrajectoryObserver.reset" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the internal state.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.GSLRandomNumberGenerator">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">GSLRandomNumberGenerator</code><a class="headerlink" href="#ecell4.core.GSLRandomNumberGenerator" title="Permalink to this definition">¶</a></dt>
<dd><p>A random number generator using the GNU Scientific Library (GSL).</p>
<p>GSLRandomNumberGenerator()</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.GSLRandomNumberGenerator.binomial">
<code class="descname">binomial</code><span class="sig-paren">(</span><em>p</em>, <em>n</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.core.GSLRandomNumberGenerator.binomial" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a random integer from the binomial distribution,
the number of successes in n independent trials with probability p.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>p</strong> : Real</p>
<blockquote>
<div><p>A probability.</p>
</div></blockquote>
<p><strong>n</strong> : Integer</p>
<blockquote>
<div><p>The number of trials.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>A random integer from a binomial distribution.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.GSLRandomNumberGenerator.gaussian">
<code class="descname">gaussian</code><span class="sig-paren">(</span><em>mean</em>, <em>sigma</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.core.GSLRandomNumberGenerator.gaussian" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a Gaussian variate with the given mean and standard deviation.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>mean</strong> : Real</p>
<blockquote>
<div><p>The mean value.</p>
</div></blockquote>
<p><strong>sigma</strong> : Real</p>
<blockquote>
<div><p>The standard deviation.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real:</p>
<blockquote class="last">
<div><p>A random number from a Gaussian distribution.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.GSLRandomNumberGenerator.seed">
<code class="descname">seed</code><span class="sig-paren">(</span><em>val=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.GSLRandomNumberGenerator.seed" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the random number seed.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>val</strong> : Integer, optional</p>
<blockquote class="last">
<div><p>A new seed. If no seed is given, reset the seed by the current time.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.GSLRandomNumberGenerator.uniform">
<code class="descname">uniform</code><span class="sig-paren">(</span><em>min</em>, <em>max</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.core.GSLRandomNumberGenerator.uniform" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a uniform random number within the given range.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>min</strong> : Real</p>
<blockquote>
<div><p>The minimum value in the range.</p>
</div></blockquote>
<p><strong>max</strong> : Real</p>
<blockquote>
<div><p>The maximum value in the range.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real:</p>
<blockquote class="last">
<div><p>A random number uniformly distributed in the range [min, max).</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.GSLRandomNumberGenerator.uniform_int">
<code class="descname">uniform_int</code><span class="sig-paren">(</span><em>min</em>, <em>max</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.core.GSLRandomNumberGenerator.uniform_int" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a uniform random number within the given range.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>min</strong> : Real</p>
<blockquote>
<div><p>The minimum value in the range.</p>
</div></blockquote>
<p><strong>max</strong> : Real</p>
<blockquote>
<div><p>The maximum value in the range.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>A random integer uniformly distributed in the range [min, max].</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Integer3">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Integer3</code><a class="headerlink" href="#ecell4.core.Integer3" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a vector consisting of three integers.</p>
<p>Integer3(Integer p1, Integer p2, Integer p3)</p>
<p class="rubric">Attributes</p>
<dl class="attribute">
<dt id="ecell4.core.Integer3.col">
<code class="descname">col</code><a class="headerlink" href="#ecell4.core.Integer3.col" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the first value.</p>
</dd></dl>

<dl class="attribute">
<dt id="ecell4.core.Integer3.layer">
<code class="descname">layer</code><a class="headerlink" href="#ecell4.core.Integer3.layer" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the third value.</p>
</dd></dl>

<dl class="attribute">
<dt id="ecell4.core.Integer3.row">
<code class="descname">row</code><a class="headerlink" href="#ecell4.core.Integer3.row" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the second value.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.MeshSurface">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">MeshSurface</code><a class="headerlink" href="#ecell4.core.MeshSurface" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a triangular mesh surface, which is
available to define structures.
The polygonal shape is given as a STL (STereoLithography) format.
This object needs VTK support.</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.MeshSurface.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.MeshSurface.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.MeshSurface.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.MeshSurface.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.MeshSurface.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.MeshSurface.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : float</p>
<blockquote class="last">
<div><p>Zero or negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Model">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Model</code><a class="headerlink" href="#ecell4.core.Model" title="Permalink to this definition">¶</a></dt>
<dd><p>A base class of a model</p>
<p>Model()</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.Model.add_reaction_rule">
<code class="descname">add_reaction_rule</code><span class="sig-paren">(</span><em>rr</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Model.add_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a new reaction rule.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>rr</strong> : ReactionRule</p>
<blockquote class="last">
<div><p>A new reaction rule.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.add_reaction_rules">
<code class="descname">add_reaction_rules</code><span class="sig-paren">(</span><em>rrs</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Model.add_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a list of new reaction rules.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>rrs</strong> : list</p>
<blockquote class="last">
<div><p>A list of new <a href="#id1"><span class="problematic" id="id2">``</span></a>ReactionRule``s.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.add_species_attribute">
<code class="descname">add_species_attribute</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Model.add_species_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a species attribute to the bottom.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote class="last">
<div><p>A new species with attributes.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.add_species_attributes">
<code class="descname">add_species_attributes</code><span class="sig-paren">(</span><em>attrs</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Model.add_species_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>Extend a list of species attributes to the bottom.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>attrs</strong> : list</p>
<blockquote class="last">
<div><p>A list of new <code class="docutils literal"><span class="pre">Species</span></code> with attributes.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.apply_species_attributes">
<code class="descname">apply_species_attributes</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Species<a class="headerlink" href="#ecell4.core.Model.apply_species_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a species with attributes.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>An original species.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Species:</p>
<blockquote class="last">
<div><p>A new species attributed by species attributes in the model.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.expand">
<code class="descname">expand</code><span class="sig-paren">(</span><em>seeds</em>, <em>max_itr=None</em>, <em>max_stoich=None</em><span class="sig-paren">)</span> &rarr; Model<a class="headerlink" href="#ecell4.core.Model.expand" title="Permalink to this definition">¶</a></dt>
<dd><p>Expand a rule-based model into a network model.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>seeds</strong> : list</p>
<blockquote>
<div><p>A list of <code class="docutils literal"><span class="pre">Species</span></code> which gives seeds.</p>
</div></blockquote>
<p><strong>max_itr</strong> : Integer</p>
<blockquote>
<div><p>A maximum number of iterations to generate new products.</p>
</div></blockquote>
<p><strong>max_stoich</strong> : Integer</p>
<blockquote>
<div><p>A maximum stoichiometry of <code class="docutils literal"><span class="pre">UnitSpecies</span></code> in a <code class="docutils literal"><span class="pre">Species</span></code>.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Model:</p>
<blockquote class="last">
<div><p>A network model.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.has_reaction_rule">
<code class="descname">has_reaction_rule</code><span class="sig-paren">(</span><em>rr</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.core.Model.has_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given reaction rule is existing or not.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.has_species_attribute">
<code class="descname">has_species_attribute</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.core.Model.has_species_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given species can be attributed or not.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.list_species">
<code class="descname">list_species</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Model.list_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of species, contained in reaction rules in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.num_reaction_rules">
<code class="descname">num_reaction_rules</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Model.num_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a number of reaction rules contained in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.query_reaction_rules">
<code class="descname">query_reaction_rules</code><span class="sig-paren">(</span><em>sp1</em>, <em>sp2=None</em><span class="sig-paren">)</span> &rarr; [ReactionRule]<a class="headerlink" href="#ecell4.core.Model.query_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Query and return a list of reaction rules, which have the given species
as their reactants.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp1</strong> : Species</p>
<blockquote>
<div><p>The first reactant</p>
</div></blockquote>
<p><strong>sp2</strong> : Species</p>
<blockquote>
<div><p>The second reactant.
This is for querying second order reaction rules.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of <a href="#id3"><span class="problematic" id="id4">``</span></a>ReactionRule``s.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.reaction_rules">
<code class="descname">reaction_rules</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Model.reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of reaction rules contained in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.remove_reaction_rule">
<code class="descname">remove_reaction_rule</code><span class="sig-paren">(</span><em>rr</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Model.remove_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove a reaction rule.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.remove_species_attribute">
<code class="descname">remove_species_attribute</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Model.remove_species_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove the species attribute.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Model.species_attributes">
<code class="descname">species_attributes</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Model.species_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of species attributes contained in the model.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.NetfreeModel">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">NetfreeModel</code><a class="headerlink" href="#ecell4.core.NetfreeModel" title="Permalink to this definition">¶</a></dt>
<dd><p>A netfree model class.</p>
<p>NetfreeModel()</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.NetfreeModel.add_parameter">
<code class="descname">add_parameter</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.add_parameter" title="Permalink to this definition">¶</a></dt>
<dd><p>This is for the tentative implementation of parameters.
This might be deprecated.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.add_parameters">
<code class="descname">add_parameters</code><span class="sig-paren">(</span><em>attrs</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.add_parameters" title="Permalink to this definition">¶</a></dt>
<dd><p>This is for the tentative implementation of parameters.
This might be deprecated.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.add_reaction_rule">
<code class="descname">add_reaction_rule</code><span class="sig-paren">(</span><em>rr</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.add_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a new reaction rule.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>rr</strong> : ReactionRule</p>
<blockquote class="last">
<div><p>A new reaction rule.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.add_reaction_rules">
<code class="descname">add_reaction_rules</code><span class="sig-paren">(</span><em>rrs</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.add_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a list of new reaction rules.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>rrs</strong> : list</p>
<blockquote class="last">
<div><p>A list of new <a href="#id5"><span class="problematic" id="id6">``</span></a>ReactionRule``s.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.add_species_attribute">
<code class="descname">add_species_attribute</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.add_species_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a species attribute to the bottom.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote class="last">
<div><p>A new species with attributes.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.add_species_attributes">
<code class="descname">add_species_attributes</code><span class="sig-paren">(</span><em>attrs</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.add_species_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>Extend a list of species attributes to the bottom.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>attrs</strong> : list</p>
<blockquote class="last">
<div><p>A list of new <code class="docutils literal"><span class="pre">Species</span></code> with attributes.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.apply_species_attributes">
<code class="descname">apply_species_attributes</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Species<a class="headerlink" href="#ecell4.core.NetfreeModel.apply_species_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a species with attributes.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>An original species.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Species:</p>
<blockquote class="last">
<div><p>A new species attributed by species attributes in the model.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.expand">
<code class="descname">expand</code><span class="sig-paren">(</span><em>seeds</em>, <em>max_itr=None</em>, <em>max_stoich=None</em><span class="sig-paren">)</span> &rarr; Model<a class="headerlink" href="#ecell4.core.NetfreeModel.expand" title="Permalink to this definition">¶</a></dt>
<dd><p>Expand a rule-based model into a network model.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>seeds</strong> : list</p>
<blockquote>
<div><p>A list of <code class="docutils literal"><span class="pre">Species</span></code> which gives seeds.</p>
</div></blockquote>
<p><strong>max_itr</strong> : Integer</p>
<blockquote>
<div><p>A maximum number of iterations to generate new products.</p>
</div></blockquote>
<p><strong>max_stoich</strong> : Integer</p>
<blockquote>
<div><p>A maximum stoichiometry of <code class="docutils literal"><span class="pre">UnitSpecies</span></code> in a <code class="docutils literal"><span class="pre">Species</span></code>.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Model:</p>
<blockquote class="last">
<div><p>A network model.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.has_reaction_rule">
<code class="descname">has_reaction_rule</code><span class="sig-paren">(</span><em>rr</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.core.NetfreeModel.has_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given reaction rule is existing or not.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.has_species_attribute">
<code class="descname">has_species_attribute</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.core.NetfreeModel.has_species_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given species can be attributed or not.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.list_species">
<code class="descname">list_species</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.list_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of species, contained in reaction rules in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.num_reaction_rules">
<code class="descname">num_reaction_rules</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.num_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a number of reaction rules contained in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.parameters">
<code class="descname">parameters</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.parameters" title="Permalink to this definition">¶</a></dt>
<dd><p>This is for the tentative implementation of parameters.
This might be deprecated.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.query_reaction_rules">
<code class="descname">query_reaction_rules</code><span class="sig-paren">(</span><em>sp1</em>, <em>sp2=None</em><span class="sig-paren">)</span> &rarr; [ReactionRule]<a class="headerlink" href="#ecell4.core.NetfreeModel.query_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Query and return a list of reaction rules, which have the given species
as their reactants.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp1</strong> : Species</p>
<blockquote>
<div><p>The first reactant</p>
</div></blockquote>
<p><strong>sp2</strong> : Species</p>
<blockquote>
<div><p>The second reactant. This is for querying second order reaction rules.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of <a href="#id7"><span class="problematic" id="id8">``</span></a>ReactionRule``s.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.reaction_rules">
<code class="descname">reaction_rules</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of reaction rules contained in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.remove_reaction_rule">
<code class="descname">remove_reaction_rule</code><span class="sig-paren">(</span><em>rr</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.remove_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove a reaction rule.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.remove_species_attribute">
<code class="descname">remove_species_attribute</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.remove_species_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove the species attribute.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetfreeModel.species_attributes">
<code class="descname">species_attributes</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetfreeModel.species_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of species attributes contained in the model.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.NetworkModel">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">NetworkModel</code><a class="headerlink" href="#ecell4.core.NetworkModel" title="Permalink to this definition">¶</a></dt>
<dd><p>A network model class.</p>
<p>NetworkModel()</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.NetworkModel.add_parameter">
<code class="descname">add_parameter</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.add_parameter" title="Permalink to this definition">¶</a></dt>
<dd><p>This is for the tentative implementation of parameters.
This might be deprecated.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.add_parameters">
<code class="descname">add_parameters</code><span class="sig-paren">(</span><em>attrs</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.add_parameters" title="Permalink to this definition">¶</a></dt>
<dd><p>This is for the tentative implementation of parameters.
This might be deprecated.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.add_reaction_rule">
<code class="descname">add_reaction_rule</code><span class="sig-paren">(</span><em>rr</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.add_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a new reaction rule.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>rr</strong> : ReactionRule</p>
<blockquote class="last">
<div><p>A new reaction rule.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.add_reaction_rules">
<code class="descname">add_reaction_rules</code><span class="sig-paren">(</span><em>rrs</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.add_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a list of new reaction rules.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>rrs</strong> : list</p>
<blockquote class="last">
<div><p>A list of new <a href="#id9"><span class="problematic" id="id10">``</span></a>ReactionRule``s.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.add_species_attribute">
<code class="descname">add_species_attribute</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.add_species_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a species attribute to the bottom.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote class="last">
<div><p>A new species with attributes.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.add_species_attributes">
<code class="descname">add_species_attributes</code><span class="sig-paren">(</span><em>attrs</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.add_species_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>Extend a list of species attributes to the bottom.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>attrs</strong> : list</p>
<blockquote class="last">
<div><p>A list of new <code class="docutils literal"><span class="pre">Species</span></code> with attributes.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.apply_species_attributes">
<code class="descname">apply_species_attributes</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Species<a class="headerlink" href="#ecell4.core.NetworkModel.apply_species_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a species with attributes.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote>
<div><p>An original species.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Species:</p>
<blockquote class="last">
<div><p>A new species attributed by species attributes in the model.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.expand">
<code class="descname">expand</code><span class="sig-paren">(</span><em>seeds</em>, <em>max_itr=None</em>, <em>max_stoich=None</em><span class="sig-paren">)</span> &rarr; Model<a class="headerlink" href="#ecell4.core.NetworkModel.expand" title="Permalink to this definition">¶</a></dt>
<dd><p>Expand a rule-based model into a network model.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>seeds</strong> : list</p>
<blockquote>
<div><p>A list of <code class="docutils literal"><span class="pre">Species</span></code> which gives seeds.</p>
</div></blockquote>
<p><strong>max_itr</strong> : Integer</p>
<blockquote>
<div><p>A maximum number of iterations to generate new products.</p>
</div></blockquote>
<p><strong>max_stoich</strong> : Integer</p>
<blockquote>
<div><p>A maximum stoichiometry of <code class="docutils literal"><span class="pre">UnitSpecies</span></code> in a <code class="docutils literal"><span class="pre">Species</span></code>.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Model:</p>
<blockquote class="last">
<div><p>A network model.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.has_reaction_rule">
<code class="descname">has_reaction_rule</code><span class="sig-paren">(</span><em>rr</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.core.NetworkModel.has_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given reaction rule is existing or not.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.has_species_attribute">
<code class="descname">has_species_attribute</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.core.NetworkModel.has_species_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given species can be attributed or not.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.list_species">
<code class="descname">list_species</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.list_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of species, contained in reaction rules in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.num_reaction_rules">
<code class="descname">num_reaction_rules</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.num_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a number of reaction rules contained in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.parameters">
<code class="descname">parameters</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.parameters" title="Permalink to this definition">¶</a></dt>
<dd><p>This is for the tentative implementation of parameters.
This might be deprecated.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.query_reaction_rules">
<code class="descname">query_reaction_rules</code><span class="sig-paren">(</span><em>sp1</em>, <em>sp2=None</em><span class="sig-paren">)</span> &rarr; [ReactionRule]<a class="headerlink" href="#ecell4.core.NetworkModel.query_reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Query and return a list of reaction rules, which have the given species
as their reactants.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp1</strong> : Species</p>
<blockquote>
<div><p>The first reactant</p>
</div></blockquote>
<p><strong>sp2</strong> : Species</p>
<blockquote>
<div><p>The second reactant. This is for querying second order reaction rules.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of <a href="#id11"><span class="problematic" id="id12">``</span></a>ReactionRule``s.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.reaction_rules">
<code class="descname">reaction_rules</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.reaction_rules" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of reaction rules contained in the model.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.remove_reaction_rule">
<code class="descname">remove_reaction_rule</code><span class="sig-paren">(</span><em>rr</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.remove_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove a reaction rule.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.remove_species_attribute">
<code class="descname">remove_species_attribute</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.remove_species_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove the species attribute.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NetworkModel.species_attributes">
<code class="descname">species_attributes</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NetworkModel.species_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of species attributes contained in the model.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.NumberObserver">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">NumberObserver</code><a class="headerlink" href="#ecell4.core.NumberObserver" title="Permalink to this definition">¶</a></dt>
<dd><p>An <code class="docutils literal"><span class="pre">Observer``class</span> <span class="pre">to</span> <span class="pre">log</span> <span class="pre">the</span> <span class="pre">number</span> <span class="pre">of</span> <span class="pre">molecules.</span>
<span class="pre">This</span> <span class="pre">``Observer</span></code> logs at the current time first, and then keeps logging
every after simulation steps.
Warning: This doesn&#8217;t work with ODESimulator.</p>
<p>NumberObserver(species)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.NumberObserver.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NumberObserver.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NumberObserver.data">
<code class="descname">data</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NumberObserver.data" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of the numbers of molecules you specified.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of lists of the number of molecules.
The size of a return value is equal to <code class="docutils literal"><span class="pre">num_steps</span></code>.
Each element of a return value is a list consisting of
time and the number of molecules specified at the construction.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NumberObserver.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NumberObserver.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the next time for logging.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NumberObserver.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NumberObserver.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NumberObserver.reset">
<code class="descname">reset</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NumberObserver.reset" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the internal state.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NumberObserver.targets">
<code class="descname">targets</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NumberObserver.targets" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of <code class="docutils literal"><span class="pre">Species</span></code>, which this <code class="docutils literal"><span class="pre">Observer</span></code> observes</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of <code class="docutils literal"><span class="pre">Species</span></code>. This is generated from arguments
you gave at the construction.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Observer">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Observer</code><a class="headerlink" href="#ecell4.core.Observer" title="Permalink to this definition">¶</a></dt>
<dd><p>A wrapper for a base class of Observers.</p>
<p>Warning: This is mainly for developers.
Do not use this for your simulation.</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.Observer.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Observer.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the next time for logging.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Observer.reset">
<code class="descname">reset</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Observer.reset" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the internal state.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Particle">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Particle</code><a class="headerlink" href="#ecell4.core.Particle" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a particle</p>
<p>Particle(Species sp, Real3 pos, Real radius, Real D)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.Particle.D">
<code class="descname">D</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Particle.D" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the diffusion coefficient.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Particle.position">
<code class="descname">position</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.core.Particle.position" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the position.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Particle.radius">
<code class="descname">radius</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Particle.radius" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the radius.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Particle.species">
<code class="descname">species</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Species<a class="headerlink" href="#ecell4.core.Particle.species" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the species.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.ParticleID">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">ParticleID</code><a class="headerlink" href="#ecell4.core.ParticleID" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing an ID of each particle</p>
<p>ParticleID(value)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.ParticleID.lot">
<code class="descname">lot</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ParticleID.lot" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the first value.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ParticleID.serial">
<code class="descname">serial</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ParticleID.serial" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the second value.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.PlanarSurface">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">PlanarSurface</code><a class="headerlink" href="#ecell4.core.PlanarSurface" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a planar surface, which is available to define
structures.</p>
<p>PlanarSurface(origin, e0, e1)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.PlanarSurface.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.PlanarSurface.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.PlanarSurface.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.PlanarSurface.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.PlanarSurface.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.PlanarSurface.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : float</p>
<blockquote class="last">
<div><p>Zero or negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.ReactionRule">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">ReactionRule</code><a class="headerlink" href="#ecell4.core.ReactionRule" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a reaction rule between <code class="docutils literal"><span class="pre">Species</span></code>.</p>
<p>ReactionRule(reactants=None, products=None, k=None)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.ReactionRule.add_product">
<code class="descname">add_product</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.add_product" title="Permalink to this definition">¶</a></dt>
<dd><p>Append a product to the end.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote class="last">
<div><p>A new product.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.add_reactant">
<code class="descname">add_reactant</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.add_reactant" title="Permalink to this definition">¶</a></dt>
<dd><p>Append a reactant to the end.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>sp</strong> : Species</p>
<blockquote class="last">
<div><p>A new reactant.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.as_string">
<code class="descname">as_string</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; str<a class="headerlink" href="#ecell4.core.ReactionRule.as_string" title="Permalink to this definition">¶</a></dt>
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
<p class="rubric">Examples</p>
<p>The string consists of a list of reactants, a list of products,
and a kinetic rate constant.</p>
<div class="highlight-python"><div class="highlight"><pre><span class="gp">&gt;&gt;&gt; </span><span class="n">rr</span> <span class="o">=</span> <span class="n">ReactionRule</span><span class="p">([</span><span class="n">Species</span><span class="p">(</span><span class="s">&quot;A&quot;</span><span class="p">),</span> <span class="n">Species</span><span class="p">(</span><span class="s">&quot;B&quot;</span><span class="p">)],</span> <span class="p">[</span><span class="n">Species</span><span class="p">(</span><span class="s">&quot;C&quot;</span><span class="p">)],</span> <span class="mf">1.0</span><span class="p">)</span>
<span class="gp">&gt;&gt;&gt; </span><span class="n">rr</span><span class="o">.</span><span class="n">as_string</span><span class="p">()</span>
<span class="go">u&#39;A+B&gt;C|1&#39;</span>
</pre></div>
</div>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.count">
<code class="descname">count</code><span class="sig-paren">(</span><em>reactants</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.core.ReactionRule.count" title="Permalink to this definition">¶</a></dt>
<dd><p>Count the number of matches for reactants.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>reactants</strong> : list</p>
<blockquote>
<div><p>A list of <code class="docutils literal"><span class="pre">Species</span></code>. The order of <code class="docutils literal"><span class="pre">reactants</span></code>
is respected.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of matches.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.ReactionRule.generate">
<code class="descname">generate</code><span class="sig-paren">(</span><em>reactants</em><span class="sig-paren">)</span> &rarr; [ReactionRule]<a class="headerlink" href="#ecell4.core.ReactionRule.generate" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate <a href="#id13"><span class="problematic" id="id14">``</span></a>ReactionRule``s from given reactants.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>reactants</strong> : list</p>
<blockquote>
<div><p>A list of <code class="docutils literal"><span class="pre">Species</span></code>. The order of <code class="docutils literal"><span class="pre">reactants</span></code> is respected.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of <code class="docutils literal"><span class="pre">ReactionRule``s.</span> <span class="pre">The</span> <span class="pre">reactants</span> <span class="pre">of</span> <span class="pre">each</span>
<span class="pre">``ReactionRule</span></code> are equal to the given <code class="docutils literal"><span class="pre">reactants</span></code>.
If the <code class="docutils literal"><span class="pre">ReactionRule</span></code> does not match the <code class="docutils literal"><span class="pre">reactants</span></code>,
return an empty list.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Examples</p>
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
<dt id="ecell4.core.ReactionRule.reactants">
<code class="descname">reactants</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.reactants" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.core.ReactionRule.set_k">
<code class="descname">set_k</code><span class="sig-paren">(</span><em>k</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.ReactionRule.set_k" title="Permalink to this definition">¶</a></dt>
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
<dt id="ecell4.core.Real3">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Real3</code><a class="headerlink" href="#ecell4.core.Real3" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a three-dimensional vector or position.</p>
<p>Real3(Real p1, Real p2, Real p3)</p>
</dd></dl>

<dl class="class">
<dt id="ecell4.core.Rod">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Rod</code><a class="headerlink" href="#ecell4.core.Rod" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a Rod shape, which is available to define
structures. The cylinder is aligned to x-axis.</p>
<p>Rod(length, radius, origin=Real3(0, 0, 0))</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.Rod.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Rod.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Rod.distance">
<code class="descname">distance</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.distance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a minimum distance from the given point to the surface.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>distance</strong> : float</p>
<blockquote class="last">
<div><p>A minimum distance from the given point.
Negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Rod.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : float</p>
<blockquote class="last">
<div><p>Zero or negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Rod.length">
<code class="descname">length</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.length" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a length of a cylinder part.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Rod.origin">
<code class="descname">origin</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.origin" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a center position of mass</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Rod.radius">
<code class="descname">radius</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.radius" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a radius of a cylinder.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Rod.shift">
<code class="descname">shift</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.shift" title="Permalink to this definition">¶</a></dt>
<dd><p>Move the center toward the given displacement</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>vec</strong> : Real3</p>
<blockquote class="last">
<div><p>A displacement.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Rod.surface">
<code class="descname">surface</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.surface" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a surface shape.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>shape</strong> : RodSurface</p>
<blockquote class="last">
<div><p>The surface shape.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.RodSurface">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">RodSurface</code><a class="headerlink" href="#ecell4.core.RodSurface" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a hollow rod surface shape, which is
available to define structures. The cylinder is aligned to x-axis.</p>
<p>RodSurface(length, radius, origin=Real3(0, 0, 0))</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.RodSurface.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.RodSurface.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.RodSurface.distance">
<code class="descname">distance</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.distance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a minimum distance from the given point to the surface.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>distance</strong> : float</p>
<blockquote class="last">
<div><p>A minimum distance from the given point.
Negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.RodSurface.inside">
<code class="descname">inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a volume shape.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>shape</strong> : Rod</p>
<blockquote class="last">
<div><p>The volume shape.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.RodSurface.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : float</p>
<blockquote class="last">
<div><p>Zero or negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.RodSurface.length">
<code class="descname">length</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.length" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a length of a cylinder part.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.RodSurface.origin">
<code class="descname">origin</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.origin" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a center position of mass</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.RodSurface.radius">
<code class="descname">radius</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.radius" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a radius of a cylinder.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.RodSurface.shift">
<code class="descname">shift</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.shift" title="Permalink to this definition">¶</a></dt>
<dd><p>Move the center toward the given displacement</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>vec</strong> : Real3</p>
<blockquote class="last">
<div><p>A displacement.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Shape">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Shape</code><a class="headerlink" href="#ecell4.core.Shape" title="Permalink to this definition">¶</a></dt>
<dd><p>A wrapper for a base class of Shapes.</p>
<p>Warning: This is mainly for developers.
Do not use this for your simulation.</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.Shape.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Shape.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Shape.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Shape.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : float</p>
<blockquote class="last">
<div><p>Zero or negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Space">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Space</code><a class="headerlink" href="#ecell4.core.Space" title="Permalink to this definition">¶</a></dt>
<dd><p>An abstract base class of all worlds. This is for developers.</p>
<p>Space()</p>
</dd></dl>

<dl class="class">
<dt id="ecell4.core.Species">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Species</code><a class="headerlink" href="#ecell4.core.Species" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a type of molecules with attributes.</p>
<p>Species(serial=None, radius=None, D=None, location=None)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.Species.add_unit">
<code class="descname">add_unit</code><span class="sig-paren">(</span><em>usp</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.add_unit" title="Permalink to this definition">¶</a></dt>
<dd><p>Append an <code class="docutils literal"><span class="pre">UnitSpecies</span></code> to the end.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>usp</strong> : UnitSpecies</p>
<blockquote class="last">
<div><p>An <code class="docutils literal"><span class="pre">UnitSpecies</span></code> to be added.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.count">
<code class="descname">count</code><span class="sig-paren">(</span><em>pttrn</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.core.Species.count" title="Permalink to this definition">¶</a></dt>
<dd><p>Count the number of matches for a pattern given as a <code class="docutils literal"><span class="pre">Species</span></code>.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pttrn</strong> : Species</p>
<blockquote>
<div><p>A pattern to be count.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of matches.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.deserialize">
<code class="descname">deserialize</code><span class="sig-paren">(</span><em>serial</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.deserialize" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the serial. All attributes will be kept.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>serial</strong> : str</p>
<blockquote class="last">
<div><p>A new serial as an unicode string.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.get_attribute">
<code class="descname">get_attribute</code><span class="sig-paren">(</span><em>name</em><span class="sig-paren">)</span> &rarr; str<a class="headerlink" href="#ecell4.core.Species.get_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Return an attribute as an unicode string.
If no corresponding attribute is found, raise an error.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>name</strong> : str</p>
<blockquote>
<div><p>The name of an attribute.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">str:</p>
<blockquote class="last">
<div><p>The value of the attribute.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.has_attribute">
<code class="descname">has_attribute</code><span class="sig-paren">(</span><em>name</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.core.Species.has_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the attribute exists or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>name</strong> : str</p>
<blockquote>
<div><p>The name of an attribute.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>True if the attribute exists, False otherwise.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.list_attributes">
<code class="descname">list_attributes</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [(str, str)]<a class="headerlink" href="#ecell4.core.Species.list_attributes" title="Permalink to this definition">¶</a></dt>
<dd><p>List all attributes.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of pairs of name and value.
<code class="docutils literal"><span class="pre">name</span></code> and <code class="docutils literal"><span class="pre">value</span></code> are given as unicode strings.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.num_units">
<code class="descname">num_units</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.core.Species.num_units" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of <code class="docutils literal"><span class="pre">UnitSpecies</span></code>.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.remove_attribute">
<code class="descname">remove_attribute</code><span class="sig-paren">(</span><em>name</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.remove_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove an attribute.
If no corresponding attribute is found, raise an error.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>name</strong> : str</p>
<blockquote class="last">
<div><p>The name of an attribute to be removed.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.serial">
<code class="descname">serial</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.serial" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the serial name as an unicode string.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.set_attribute">
<code class="descname">set_attribute</code><span class="sig-paren">(</span><em>name</em>, <em>value</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Species.set_attribute" title="Permalink to this definition">¶</a></dt>
<dd><p>Set an attribute.
If existing already, the attribute will be overwritten.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>name</strong> : str</p>
<blockquote>
<div><p>The name of an attribute.</p>
</div></blockquote>
<p><strong>value</strong> : str</p>
<blockquote class="last">
<div><p>The value of an attribute.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Species.units">
<code class="descname">units</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; [UnitSpecies]<a class="headerlink" href="#ecell4.core.Species.units" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of all <code class="docutils literal"><span class="pre">UnitSpecies</span></code> contained.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Sphere">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Sphere</code><a class="headerlink" href="#ecell4.core.Sphere" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a sphere shape, which is available to define
structures.</p>
<p>Sphere(center, radius)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.Sphere.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Sphere.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Sphere.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Sphere.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Sphere.distance">
<code class="descname">distance</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Sphere.distance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a minimum distance from the given point to the surface.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>distance</strong> : float</p>
<blockquote class="last">
<div><p>A minimum distance from the given point.
Negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Sphere.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Sphere.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : float</p>
<blockquote class="last">
<div><p>Zero or negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Sphere.surface">
<code class="descname">surface</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Sphere.surface" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a surface shape.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>shape</strong> : SphericalSurface</p>
<blockquote class="last">
<div><p>The surface shape.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.SphericalSurface">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">SphericalSurface</code><a class="headerlink" href="#ecell4.core.SphericalSurface" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a hollow spherical surface, which is
available to define structures.</p>
<p>SphericalSurface(center, radius)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.SphericalSurface.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.SphericalSurface.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.SphericalSurface.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.SphericalSurface.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.SphericalSurface.distance">
<code class="descname">distance</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.SphericalSurface.distance" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a minimum distance from the given point to the surface.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>distance</strong> : float</p>
<blockquote class="last">
<div><p>A minimum distance from the given point.
Negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.SphericalSurface.inside">
<code class="descname">inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.SphericalSurface.inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a volume shape.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>shape</strong> : Sphere</p>
<blockquote class="last">
<div><p>The volume shape.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.SphericalSurface.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.SphericalSurface.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pos</strong> : Real3</p>
<blockquote>
<div><p>A position.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first"><strong>value</strong> : float</p>
<blockquote class="last">
<div><p>Zero or negative if the given point is inside.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.TimeoutObserver">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">TimeoutObserver</code><a class="headerlink" href="#ecell4.core.TimeoutObserver" title="Permalink to this definition">¶</a></dt>
<dd><p>An <a href="#id15"><span class="problematic" id="id16">``</span></a>Observer``class to stop simulation at the given calculation time.</p>
<p>TimeoutObserver(interval)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.TimeoutObserver.accumulation">
<code class="descname">accumulation</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimeoutObserver.accumulation" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the accumulation time.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.TimeoutObserver.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimeoutObserver.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.TimeoutObserver.duration">
<code class="descname">duration</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimeoutObserver.duration" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the last time to be called.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.TimeoutObserver.interval">
<code class="descname">interval</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimeoutObserver.interval" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the timeout in seconds.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.TimeoutObserver.reset">
<code class="descname">reset</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimeoutObserver.reset" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the internal state.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.TimingNumberObserver">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">TimingNumberObserver</code><a class="headerlink" href="#ecell4.core.TimingNumberObserver" title="Permalink to this definition">¶</a></dt>
<dd><p>An <a href="#id17"><span class="problematic" id="id18">``</span></a>Observer``class to log the number of molecules just at the time
you assigned.</p>
<p>TimingNumberObserver(t, species)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.TimingNumberObserver.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimingNumberObserver.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.TimingNumberObserver.data">
<code class="descname">data</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimingNumberObserver.data" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of the numbers of molecules you specified.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of lists of the number of molecules.
The size of a return value is equal to <code class="docutils literal"><span class="pre">num_steps</span></code>.
Each element of a return value is a list consisting of
time and the number of molecules specified at the construction.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.TimingNumberObserver.next_time">
<code class="descname">next_time</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimingNumberObserver.next_time" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the next time for logging.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.TimingNumberObserver.num_steps">
<code class="descname">num_steps</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimingNumberObserver.num_steps" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of steps.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.TimingNumberObserver.reset">
<code class="descname">reset</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimingNumberObserver.reset" title="Permalink to this definition">¶</a></dt>
<dd><p>Reset the internal state.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.TimingNumberObserver.targets">
<code class="descname">targets</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimingNumberObserver.targets" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of <code class="docutils literal"><span class="pre">Species</span></code>, which this <code class="docutils literal"><span class="pre">Observer</span></code> observes</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of <code class="docutils literal"><span class="pre">Species</span></code>. This is generated from arguments
you gave at the construction.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.UnitSpecies">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">UnitSpecies</code><a class="headerlink" href="#ecell4.core.UnitSpecies" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing an unit of species.</p>
<p>UnitSpecies(name=None)</p>
<div class="admonition seealso">
<p class="first admonition-title">See also</p>
<p class="last"><a class="reference internal" href="#ecell4.core.Species" title="ecell4.core.Species"><code class="xref py py-obj docutils literal"><span class="pre">Species</span></code></a></p>
</div>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.UnitSpecies.add_site">
<code class="descname">add_site</code><span class="sig-paren">(</span><em>name</em>, <em>state</em>, <em>bond</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.UnitSpecies.add_site" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a new site.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>name</strong> : str</p>
<blockquote>
<div><p>A name of the site</p>
</div></blockquote>
<p><strong>state</strong> : str</p>
<blockquote>
<div><p>A state name of the site</p>
</div></blockquote>
<p><strong>bond</strong> : str</p>
<blockquote class="last">
<div><p>A bond of the site.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.UnitSpecies.deserialize">
<code class="descname">deserialize</code><span class="sig-paren">(</span><em>serial</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.UnitSpecies.deserialize" title="Permalink to this definition">¶</a></dt>
<dd><p>Deserialize the given serial, and load it.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>serial</strong> : str</p>
<blockquote class="last">
<div><p>A serial</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.UnitSpecies.name">
<code class="descname">name</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.UnitSpecies.name" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a name.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.UnitSpecies.serial">
<code class="descname">serial</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.UnitSpecies.serial" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the serial, which consists of a name and sites.</p>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Voxel">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Voxel</code><a class="headerlink" href="#ecell4.core.Voxel" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a voxel in LatticeSpace.</p>
<p>Voxel(Species sp, Integer coord, Real radius, Real D, loc=None)</p>
<p class="rubric">Methods</p>
<dl class="method">
<dt id="ecell4.core.Voxel.D">
<code class="descname">D</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Voxel.D" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the diffusion coefficient.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Voxel.coordinate">
<code class="descname">coordinate</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Voxel.coordinate" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the coordinate.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Voxel.loc">
<code class="descname">loc</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; str<a class="headerlink" href="#ecell4.core.Voxel.loc" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the location information as a string.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Voxel.radius">
<code class="descname">radius</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Voxel.radius" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the radius.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Voxel.species">
<code class="descname">species</code><span class="sig-paren">(</span><span class="sig-paren">)</span> &rarr; Species<a class="headerlink" href="#ecell4.core.Voxel.species" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the species.</p>
</dd></dl>

</dd></dl>

<dl class="function">
<dt id="ecell4.core.cbrt">
<code class="descclassname">ecell4.core.</code><code class="descname">cbrt</code><span class="sig-paren">(</span><em>x</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.core.cbrt" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a cubic root of the given value.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.count_rrmatches">
<code class="descclassname">ecell4.core.</code><code class="descname">count_rrmatches</code><span class="sig-paren">(</span><em>pttrn</em>, <em>reactants</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.core.count_rrmatches" title="Permalink to this definition">¶</a></dt>
<dd><p>Count the number of matches for a pattern given as a <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pttrn</strong> : ReactionRule</p>
<blockquote>
<div><p>A pattern.</p>
</div></blockquote>
<p><strong>reactants</strong> : list</p>
<blockquote>
<div><p>A list of reactants, <code class="docutils literal"><span class="pre">Species</span></code>. The order of reactants is respected.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of matches.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.count_spmatches">
<code class="descclassname">ecell4.core.</code><code class="descname">count_spmatches</code><span class="sig-paren">(</span><em>pttrn</em>, <em>sp</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.core.count_spmatches" title="Permalink to this definition">¶</a></dt>
<dd><p>Count the number of matches for a pattern given as a <code class="docutils literal"><span class="pre">Species</span></code>.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pttrn</strong> : Species</p>
<blockquote>
<div><p>A pattern.</p>
</div></blockquote>
<p><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A target.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer:</p>
<blockquote class="last">
<div><p>The number of matches.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Notes</p>
<p>Rather use <code class="docutils literal"><span class="pre">Species.count</span></code>.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.create_binding_reaction_rule">
<code class="descclassname">ecell4.core.</code><code class="descname">create_binding_reaction_rule</code><span class="sig-paren">(</span><em>reactant1</em>, <em>reactant2</em>, <em>product1</em>, <em>k</em><span class="sig-paren">)</span> &rarr; ReactionRule<a class="headerlink" href="#ecell4.core.create_binding_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Create a binding <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>reactant1</strong> : Species</p>
<blockquote>
<div><p>One of two reactants.</p>
</div></blockquote>
<p><strong>reactant2</strong> : Species</p>
<blockquote>
<div><p>One of two reactants.</p>
</div></blockquote>
<p><strong>product1</strong> : Species</p>
<blockquote>
<div><p>A product.</p>
</div></blockquote>
<p><strong>k</strong> : float</p>
<blockquote class="last">
<div><p>A kinetic parameter.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Notes</p>
<p>This is equivalent to <code class="docutils literal"><span class="pre">ReactionRule([reactant1,</span> <span class="pre">reactant2],</span> <span class="pre">[product1],</span> <span class="pre">k)</span></code>.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.create_degradation_reaction_rule">
<code class="descclassname">ecell4.core.</code><code class="descname">create_degradation_reaction_rule</code><span class="sig-paren">(</span><em>reactant1</em>, <em>k</em><span class="sig-paren">)</span> &rarr; ReactionRule<a class="headerlink" href="#ecell4.core.create_degradation_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Create a degradation <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>reactant1</strong> : Species</p>
<blockquote>
<div><p>A reactant to be degradated.</p>
</div></blockquote>
<p><strong>k</strong> : float</p>
<blockquote class="last">
<div><p>A kinetic parameter.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Notes</p>
<p>This is equivalent to <code class="docutils literal"><span class="pre">ReactionRule([reactant1],</span> <span class="pre">[],</span> <span class="pre">k)</span></code>.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.create_synthesis_reaction_rule">
<code class="descclassname">ecell4.core.</code><code class="descname">create_synthesis_reaction_rule</code><span class="sig-paren">(</span><em>product1</em>, <em>k</em><span class="sig-paren">)</span> &rarr; ReactionRule<a class="headerlink" href="#ecell4.core.create_synthesis_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Create a synthesis <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>product1</strong> : Species</p>
<blockquote>
<div><p>A product to be synthesized.</p>
</div></blockquote>
<p><strong>k</strong> : float</p>
<blockquote class="last">
<div><p>A kinetic parameter.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Notes</p>
<p>This is equivalent to <code class="docutils literal"><span class="pre">ReactionRule([],</span> <span class="pre">[product1],</span> <span class="pre">k)</span></code>.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.create_unbinding_reaction_rule">
<code class="descclassname">ecell4.core.</code><code class="descname">create_unbinding_reaction_rule</code><span class="sig-paren">(</span><em>reactant1</em>, <em>product1</em>, <em>product2</em>, <em>k</em><span class="sig-paren">)</span> &rarr; ReactionRule<a class="headerlink" href="#ecell4.core.create_unbinding_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>Create an unbinding <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>reactant1</strong> : Species</p>
<blockquote>
<div><p>A reactant.</p>
</div></blockquote>
<p><strong>product1</strong> : Species</p>
<blockquote>
<div><p>One of two products.</p>
</div></blockquote>
<p><strong>product2</strong> : Species</p>
<blockquote>
<div><p>One of two products.</p>
</div></blockquote>
<p><strong>k</strong> : float</p>
<blockquote class="last">
<div><p>A kinetic parameter.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Notes</p>
<p>This is equivalent to <code class="docutils literal"><span class="pre">ReactionRule([reactant1],</span> <span class="pre">[product1,</span> <span class="pre">product2],</span> <span class="pre">k)</span></code>.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.create_unimolecular_reaction_rule">
<code class="descclassname">ecell4.core.</code><code class="descname">create_unimolecular_reaction_rule</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.create_unimolecular_reaction_rule" title="Permalink to this definition">¶</a></dt>
<dd><p>create_synthesis_reaction_rule(reactant1, product1, k) -&gt; ReactionRule</p>
<p>Create an unimolecular <code class="docutils literal"><span class="pre">ReactionRule</span></code>.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>reactant1</strong> : Species</p>
<blockquote>
<div><p>A reactant to be modified.</p>
</div></blockquote>
<p><strong>product1</strong> : Species</p>
<blockquote>
<div><p>A product.</p>
</div></blockquote>
<p><strong>k</strong> : float</p>
<blockquote class="last">
<div><p>A kinetic parameter.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Notes</p>
<p>This is equivalent to <code class="docutils literal"><span class="pre">ReactionRule([reactant1],</span> <span class="pre">[product1],</span> <span class="pre">k)</span></code>.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.cross_product">
<code class="descclassname">ecell4.core.</code><code class="descname">cross_product</code><span class="sig-paren">(</span><em>p1</em>, <em>p2</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.core.cross_product" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a cross product between two vectors</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.dot_product">
<code class="descclassname">ecell4.core.</code><code class="descname">dot_product</code><span class="sig-paren">(</span><em>p1</em>, <em>p2</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.core.dot_product" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dot product between two vectors</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.format_species">
<code class="descclassname">ecell4.core.</code><code class="descname">format_species</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Species<a class="headerlink" href="#ecell4.core.format_species" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a species uniquely reformatted.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.integer3_abs">
<code class="descclassname">ecell4.core.</code><code class="descname">integer3_abs</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.integer3_abs" title="Permalink to this definition">¶</a></dt>
<dd><p>real3_abs(p1) -&gt; Integer3</p>
<p>Return an absolute vector of the given vector.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>p1</strong> : Integer3</p>
<blockquote>
<div><p>A vector.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer3:</p>
<blockquote class="last">
<div><p>The absolute vector, which consists of absolute value
of the given vector.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<div class="admonition seealso">
<p class="first admonition-title">See also</p>
<p class="last"><a class="reference internal" href="#ecell4.core.length" title="ecell4.core.length"><code class="xref py py-obj docutils literal"><span class="pre">length</span></code></a></p>
</div>
<p class="rubric">Notes</p>
<p>This is NOT for taking the norm of a vector.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.integer3_add">
<code class="descclassname">ecell4.core.</code><code class="descname">integer3_add</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.integer3_add" title="Permalink to this definition">¶</a></dt>
<dd><p>add(p1, p2) -&gt; Integer3</p>
<p>Add two <a href="#id19"><span class="problematic" id="id20">``</span></a>Integer3``s, and returns the sum.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>p1</strong> : Integer3</p>
<blockquote>
<div><p>The first vector.</p>
</div></blockquote>
<p><strong>p2</strong> : Integer3</p>
<blockquote>
<div><p>The second vector.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer3:</p>
<blockquote class="last">
<div><p>The sum of two vectors, <code class="docutils literal"><span class="pre">p1</span> <span class="pre">+</span> <span class="pre">p2</span></code>.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.integer3_subtract">
<code class="descclassname">ecell4.core.</code><code class="descname">integer3_subtract</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.integer3_subtract" title="Permalink to this definition">¶</a></dt>
<dd><p>subtract(p1, p2) -&gt; Integer3</p>
<p>Subtract p2 from p1.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>p1</strong> : Integer3</p>
<blockquote>
<div><p>The left-hand-side vector.</p>
</div></blockquote>
<p><strong>p2</strong> : Integer3</p>
<blockquote>
<div><p>The right-hand-side vector.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Integer3:</p>
<blockquote class="last">
<div><p>Its difference, <code class="docutils literal"><span class="pre">p1</span> <span class="pre">-</span> <span class="pre">p2</span></code>.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.length">
<code class="descclassname">ecell4.core.</code><code class="descname">length</code><span class="sig-paren">(</span><em>p1</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.core.length" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a Euclidean norm of the given vector.
This is almost equivalent to call <code class="docutils literal"><span class="pre">sqrt(length_sq(p1))</span></code></p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.length_sq">
<code class="descclassname">ecell4.core.</code><code class="descname">length_sq</code><span class="sig-paren">(</span><em>p1</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.core.length_sq" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a square of a Euclidean norm of the given vector.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.load_version_information">
<code class="descclassname">ecell4.core.</code><code class="descname">load_version_information</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.load_version_information" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a version information of HDF5 as a string.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.real3_abs">
<code class="descclassname">ecell4.core.</code><code class="descname">real3_abs</code><span class="sig-paren">(</span><em>p1</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.core.real3_abs" title="Permalink to this definition">¶</a></dt>
<dd><p>Return an absolute vector of the given vector.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>p1</strong> : Real3</p>
<blockquote>
<div><p>A vector.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real3:</p>
<blockquote class="last">
<div><p>The absolute vector, which consists of absolute value of the given vector.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<div class="admonition seealso">
<p class="first admonition-title">See also</p>
<p class="last"><a class="reference internal" href="#ecell4.core.length" title="ecell4.core.length"><code class="xref py py-obj docutils literal"><span class="pre">length</span></code></a></p>
</div>
<p class="rubric">Notes</p>
<p>This is NOT for taking the norm of a vector.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.real3_add">
<code class="descclassname">ecell4.core.</code><code class="descname">real3_add</code><span class="sig-paren">(</span><em>p1</em>, <em>p2</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.core.real3_add" title="Permalink to this definition">¶</a></dt>
<dd><p>Add two <a href="#id21"><span class="problematic" id="id22">``</span></a>Real3``s, and returns the sum.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>p1</strong> : Real3</p>
<blockquote>
<div><p>The first vector.</p>
</div></blockquote>
<p><strong>p2</strong> : Real3</p>
<blockquote>
<div><p>The second vector.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real3:</p>
<blockquote class="last">
<div><p>The sum of two vectors, <code class="docutils literal"><span class="pre">p1</span> <span class="pre">+</span> <span class="pre">p2</span></code>.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.real3_divide">
<code class="descclassname">ecell4.core.</code><code class="descname">real3_divide</code><span class="sig-paren">(</span><em>p1</em>, <em>p2</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.core.real3_divide" title="Permalink to this definition">¶</a></dt>
<dd><p>Divide p1 by p2.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>p1</strong> : Real3</p>
<blockquote>
<div><p>The numerator.</p>
</div></blockquote>
<p><strong>p2</strong> : Real</p>
<blockquote>
<div><p>The denominator.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real3:</p>
<blockquote class="last">
<div><p>The divided vector, <code class="docutils literal"><span class="pre">p1</span> <span class="pre">/</span> <span class="pre">p2</span></code>.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.real3_multiply">
<code class="descclassname">ecell4.core.</code><code class="descname">real3_multiply</code><span class="sig-paren">(</span><em>p1</em>, <em>p2</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.core.real3_multiply" title="Permalink to this definition">¶</a></dt>
<dd><p>Multiply p1 by p2.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>p1</strong> : Real3</p>
<blockquote>
<div><p>A vector.</p>
</div></blockquote>
<p><strong>p2</strong> : Real</p>
<blockquote>
<div><p>A factor.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real3:</p>
<blockquote class="last">
<div><p>The multipled vector, <code class="docutils literal"><span class="pre">p1</span> <span class="pre">*</span> <span class="pre">p2</span></code>.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.real3_subtract">
<code class="descclassname">ecell4.core.</code><code class="descname">real3_subtract</code><span class="sig-paren">(</span><em>p1</em>, <em>p2</em><span class="sig-paren">)</span> &rarr; Real3<a class="headerlink" href="#ecell4.core.real3_subtract" title="Permalink to this definition">¶</a></dt>
<dd><p>Subtract p2 from p1.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>p1</strong> : Real3</p>
<blockquote>
<div><p>The left-hand-side vector.</p>
</div></blockquote>
<p><strong>p2</strong> : Real3</p>
<blockquote>
<div><p>The right-hand-side vector.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">Real3:</p>
<blockquote class="last">
<div><p>Its difference, <code class="docutils literal"><span class="pre">p1</span> <span class="pre">-</span> <span class="pre">p2</span></code>.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.rrgenerate">
<code class="descclassname">ecell4.core.</code><code class="descname">rrgenerate</code><span class="sig-paren">(</span><em>pttrn</em>, <em>reactants</em><span class="sig-paren">)</span> &rarr; [Species]<a class="headerlink" href="#ecell4.core.rrgenerate" title="Permalink to this definition">¶</a></dt>
<dd><p>Generate a list of products from the given list of reactants.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pttrn</strong> : ReactionRule</p>
<blockquote>
<div><p>A pattern.</p>
</div></blockquote>
<p><strong>reactants</strong> : list</p>
<blockquote>
<div><p>A list of <code class="docutils literal"><span class="pre">Species</span></code>. The order of <code class="docutils literal"><span class="pre">reactants</span></code> is respected.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">list:</p>
<blockquote class="last">
<div><p>A list of products. The size of the list is equal to the number of matches.
Each element of the list is a list of <code class="docutils literal"><span class="pre">Species</span></code>.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
<p class="rubric">Notes</p>
<p>Rather use <code class="docutils literal"><span class="pre">ReactionRule.generate</span></code>.</p>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.rrmatch">
<code class="descclassname">ecell4.core.</code><code class="descname">rrmatch</code><span class="sig-paren">(</span><em>pttrn</em>, <em>reactants</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.core.rrmatch" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if a pattern matches the reactants or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pttrn</strong> : ReactionRule</p>
<blockquote>
<div><p>A pattern.</p>
</div></blockquote>
<p><strong>reactants</strong> : list</p>
<blockquote>
<div><p>A list of reactants, <code class="docutils literal"><span class="pre">Species</span></code>. The order of reactants is respected.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>True if <code class="docutils literal"><span class="pre">pttrn</span></code> matches <code class="docutils literal"><span class="pre">reactants</span></code> at least one time,
False otherwise.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.spmatch">
<code class="descclassname">ecell4.core.</code><code class="descname">spmatch</code><span class="sig-paren">(</span><em>pttrn</em>, <em>sp</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.core.spmatch" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if a pattern matches the target <code class="docutils literal"><span class="pre">Species</span></code> or not.</p>
<table class="docutils field-list" frame="void" rules="none">
<col class="field-name" />
<col class="field-body" />
<tbody valign="top">
<tr class="field-odd field"><th class="field-name">Parameters:</th><td class="field-body"><p class="first"><strong>pttrn</strong> : Species</p>
<blockquote>
<div><p>A pattern.</p>
</div></blockquote>
<p><strong>sp</strong> : Species</p>
<blockquote>
<div><p>A target.</p>
</div></blockquote>
</td>
</tr>
<tr class="field-even field"><th class="field-name">Returns:</th><td class="field-body"><p class="first">bool:</p>
<blockquote class="last">
<div><p>True if <code class="docutils literal"><span class="pre">pttrn</span></code> matches <code class="docutils literal"><span class="pre">sp</span></code> at least one time, False otherwise.</p>
</div></blockquote>
</td>
</tr>
</tbody>
</table>
</dd></dl>

<dl class="function">
<dt id="ecell4.core.unique_serial">
<code class="descclassname">ecell4.core.</code><code class="descname">unique_serial</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; str<a class="headerlink" href="#ecell4.core.unique_serial" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a serial of a species uniquely reformatted.
This is equivalent to call <code class="docutils literal"><span class="pre">format_species(sp).serial()</span></code></p>
</dd></dl>
