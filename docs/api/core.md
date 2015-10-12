  <span class="target" id="module-ecell4.core"></span><dl class="class">
<dt id="ecell4.core.AABB">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">AABB</code><a class="headerlink" href="#ecell4.core.AABB" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing an axis aligned bounding box (AABB),
which is available to define structures.</p>
<p>AABB(lower, upper)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.AABB.as_base" title="ecell4.core.AABB.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.AABB.dimension" title="ecell4.core.AABB.dimension"><code class="xref py py-obj docutils literal"><span class="pre">dimension</span></code></a></td>
<td>Return a dimension of this shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.AABB.distance" title="ecell4.core.AABB.distance"><code class="xref py py-obj docutils literal"><span class="pre">distance</span></code></a></td>
<td>Return a minimum distance from the given point to the surface.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.AABB.is_inside" title="ecell4.core.AABB.is_inside"><code class="xref py py-obj docutils literal"><span class="pre">is_inside</span></code></a></td>
<td>Return if the given point is inside or not.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.AABB.lower" title="ecell4.core.AABB.lower"><code class="xref py py-obj docutils literal"><span class="pre">lower</span></code></a></td>
<td>Return a vertex suggesting the lower bounds.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.AABB.upper" title="ecell4.core.AABB.upper"><code class="xref py py-obj docutils literal"><span class="pre">upper</span></code></a></td>
<td>Return a vertex suggesting the upper bounds.</td>
</tr>
</tbody>
</table>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>distance (float): A minimum distance from the given point.</dt>
<dd>Negative if the given point is inside.</dd>
</dl>
</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.AABB.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.AABB.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd>value (float): Zero or negative if the given point is inside.</dd>
</dl>
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
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Cylinder.as_base" title="ecell4.core.Cylinder.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Cylinder.dimension" title="ecell4.core.Cylinder.dimension"><code class="xref py py-obj docutils literal"><span class="pre">dimension</span></code></a></td>
<td>Return a dimension of this shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Cylinder.distance" title="ecell4.core.Cylinder.distance"><code class="xref py py-obj docutils literal"><span class="pre">distance</span></code></a></td>
<td>Return a minimum distance from the given point to the surface.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Cylinder.is_inside" title="ecell4.core.Cylinder.is_inside"><code class="xref py py-obj docutils literal"><span class="pre">is_inside</span></code></a></td>
<td>Return if the given point is inside or not.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Cylinder.surface" title="ecell4.core.Cylinder.surface"><code class="xref py py-obj docutils literal"><span class="pre">surface</span></code></a></td>
<td>Create and return a surface shape.</td>
</tr>
</tbody>
</table>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>distance (float): A minimum distance from the given point.</dt>
<dd>Negative if the given point is inside.</dd>
</dl>
</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Cylinder.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Cylinder.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd>value (float): Zero or negative if the given point is inside.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Cylinder.surface">
<code class="descname">surface</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Cylinder.surface" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a surface shape.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd>shape (CylindricalSurface): The surface shape.</dd>
</dl>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.CylindricalSurface">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">CylindricalSurface</code><a class="headerlink" href="#ecell4.core.CylindricalSurface" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a hollow cylindrical surface, which is
available to define structures.</p>
<p>CylindricalSurface(center, radius, axis, half_height)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.CylindricalSurface.as_base" title="ecell4.core.CylindricalSurface.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.CylindricalSurface.dimension" title="ecell4.core.CylindricalSurface.dimension"><code class="xref py py-obj docutils literal"><span class="pre">dimension</span></code></a></td>
<td>Return a dimension of this shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.CylindricalSurface.distance" title="ecell4.core.CylindricalSurface.distance"><code class="xref py py-obj docutils literal"><span class="pre">distance</span></code></a></td>
<td>Return a minimum distance from the given point to the surface.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.CylindricalSurface.inside" title="ecell4.core.CylindricalSurface.inside"><code class="xref py py-obj docutils literal"><span class="pre">inside</span></code></a></td>
<td>Create and return a volume shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.CylindricalSurface.is_inside" title="ecell4.core.CylindricalSurface.is_inside"><code class="xref py py-obj docutils literal"><span class="pre">is_inside</span></code></a></td>
<td>Return if the given point is inside or not.</td>
</tr>
</tbody>
</table>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>distance (float): A minimum distance from the given point.</dt>
<dd>Negative if the given point is inside.</dd>
</dl>
</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.CylindricalSurface.inside">
<code class="descname">inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.CylindricalSurface.inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a volume shape.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd>shape (Cylinder): The volume shape.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.CylindricalSurface.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.CylindricalSurface.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd>value (float): Zero or negative if the given point is inside.</dd>
</dl>
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
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalCSVObserver.as_base" title="ecell4.core.FixedIntervalCSVObserver.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.FixedIntervalCSVObserver.filename" title="ecell4.core.FixedIntervalCSVObserver.filename"><code class="xref py py-obj docutils literal"><span class="pre">filename</span></code></a></td>
<td>Return a file name to be saved at the next time</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalCSVObserver.log" title="ecell4.core.FixedIntervalCSVObserver.log"><code class="xref py py-obj docutils literal"><span class="pre">log</span></code></a></td>
<td>Force to log the given <code class="docutils literal"><span class="pre">World</span></code> to a file.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.FixedIntervalCSVObserver.next_time" title="ecell4.core.FixedIntervalCSVObserver.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the next time for logging.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalCSVObserver.num_steps" title="ecell4.core.FixedIntervalCSVObserver.num_steps"><code class="xref py py-obj docutils literal"><span class="pre">num_steps</span></code></a></td>
<td>Return the number of steps.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.FixedIntervalCSVObserver.reset" title="ecell4.core.FixedIntervalCSVObserver.reset"><code class="xref py py-obj docutils literal"><span class="pre">reset</span></code></a></td>
<td>Reset the internal state.</td>
</tr>
</tbody>
</table>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>w (Space): A <code class="docutils literal"><span class="pre">Space</span></code> (<code class="docutils literal"><span class="pre">World</span></code>) to be logged.</dd>
<dt>Example:</dt>
<dd><p class="first">This is an easy way to save a <code class="docutils literal"><span class="pre">World</span></code> in CSV format without
running a simulation.</p>
<div class="last highlight-python"><div class="highlight"><pre><span class="gp">&gt;&gt;&gt; </span><span class="n">w</span> <span class="o">=</span> <span class="n">lattice</span><span class="o">.</span><span class="n">LatticeWorld</span><span class="p">(</span><span class="n">Real3</span><span class="p">(</span><span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">,</span> <span class="mi">1</span><span class="p">),</span> <span class="mf">0.005</span><span class="p">)</span>
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
</dd>
</dl>
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
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalHDF5Observer.as_base" title="ecell4.core.FixedIntervalHDF5Observer.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.FixedIntervalHDF5Observer.filename" title="ecell4.core.FixedIntervalHDF5Observer.filename"><code class="xref py py-obj docutils literal"><span class="pre">filename</span></code></a></td>
<td>Return a file name to be saved at the next time</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalHDF5Observer.next_time" title="ecell4.core.FixedIntervalHDF5Observer.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the next time for logging.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.FixedIntervalHDF5Observer.num_steps" title="ecell4.core.FixedIntervalHDF5Observer.num_steps"><code class="xref py py-obj docutils literal"><span class="pre">num_steps</span></code></a></td>
<td>Return the number of steps.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalHDF5Observer.reset" title="ecell4.core.FixedIntervalHDF5Observer.reset"><code class="xref py py-obj docutils literal"><span class="pre">reset</span></code></a></td>
<td>Reset the internal state.</td>
</tr>
</tbody>
</table>
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
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalNumberObserver.as_base" title="ecell4.core.FixedIntervalNumberObserver.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.FixedIntervalNumberObserver.data" title="ecell4.core.FixedIntervalNumberObserver.data"><code class="xref py py-obj docutils literal"><span class="pre">data</span></code></a></td>
<td>Return a list of the number of molecules you specified.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalNumberObserver.next_time" title="ecell4.core.FixedIntervalNumberObserver.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the next time for logging.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.FixedIntervalNumberObserver.num_steps" title="ecell4.core.FixedIntervalNumberObserver.num_steps"><code class="xref py py-obj docutils literal"><span class="pre">num_steps</span></code></a></td>
<td>Return the number of steps.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalNumberObserver.reset" title="ecell4.core.FixedIntervalNumberObserver.reset"><code class="xref py py-obj docutils literal"><span class="pre">reset</span></code></a></td>
<td>Reset the internal state.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.FixedIntervalNumberObserver.targets" title="ecell4.core.FixedIntervalNumberObserver.targets"><code class="xref py py-obj docutils literal"><span class="pre">targets</span></code></a></td>
<td>Return a list of <code class="docutils literal"><span class="pre">Species</span></code>, which this <code class="docutils literal"><span class="pre">Observer</span></code> observes</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.core.FixedIntervalNumberObserver.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalNumberObserver.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalNumberObserver.data">
<code class="descname">data</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalNumberObserver.data" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of the number of molecules you specified.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of lists of the numbers of molecules.</dt>
<dd>The size of a return value is equal to <code class="docutils literal"><span class="pre">num_steps</span></code>.
Each element of a return value is a list consisting of
time and the number of molecules specified at the construction.</dd>
</dl>
</dd>
</dl>
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
<dl class="docutils">
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of <code class="docutils literal"><span class="pre">Species</span></code>. This is generated from arguments</dt>
<dd>you gave at the construction.</dd>
</dl>
</dd>
</dl>
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
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalTrajectoryObserver.as_base" title="ecell4.core.FixedIntervalTrajectoryObserver.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.FixedIntervalTrajectoryObserver.data" title="ecell4.core.FixedIntervalTrajectoryObserver.data"><code class="xref py py-obj docutils literal"><span class="pre">data</span></code></a></td>
<td>Return a list of trajectories for each particles.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalTrajectoryObserver.next_time" title="ecell4.core.FixedIntervalTrajectoryObserver.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the next time for logging.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.FixedIntervalTrajectoryObserver.num_steps" title="ecell4.core.FixedIntervalTrajectoryObserver.num_steps"><code class="xref py py-obj docutils literal"><span class="pre">num_steps</span></code></a></td>
<td>Return the number of steps.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.FixedIntervalTrajectoryObserver.reset" title="ecell4.core.FixedIntervalTrajectoryObserver.reset"><code class="xref py py-obj docutils literal"><span class="pre">reset</span></code></a></td>
<td>Reset the internal state.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.core.FixedIntervalTrajectoryObserver.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalTrajectoryObserver.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.FixedIntervalTrajectoryObserver.data">
<code class="descname">data</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.FixedIntervalTrajectoryObserver.data" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of trajectories for each particles.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of lists of <code class="docutils literal"><span class="pre">Real3</span></code>. An element of a return value</dt>
<dd>is corresponding the trajectory of each particle. Thus, the size
of a return value is the same with that of <code class="docutils literal"><span class="pre">pids</span></code> you gave
at the construction.
If a particle corresponding to the given <code class="docutils literal"><span class="pre">ParticleID</span></code> is missing,
i.e. for a reaction, this <code class="docutils literal"><span class="pre">Observer</span></code> just skips to log the
position. Therefore, lengths of the trajectories can be diverse.</dd>
</dl>
</dd>
</dl>
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
<dt id="ecell4.core.MeshSurface">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">MeshSurface</code><a class="headerlink" href="#ecell4.core.MeshSurface" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a triangular mesh surface, which is
available to define structures.
The polygonal shape is given as a STL (STereoLithography) format.
This object needs VTK support.</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.MeshSurface.as_base" title="ecell4.core.MeshSurface.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.MeshSurface.dimension" title="ecell4.core.MeshSurface.dimension"><code class="xref py py-obj docutils literal"><span class="pre">dimension</span></code></a></td>
<td>Return a dimension of this shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.MeshSurface.is_inside" title="ecell4.core.MeshSurface.is_inside"><code class="xref py py-obj docutils literal"><span class="pre">is_inside</span></code></a></td>
<td>Return if the given point is inside or not.</td>
</tr>
</tbody>
</table>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd>value (float): Zero or negative if the given point is inside.</dd>
</dl>
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
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.NumberObserver.as_base" title="ecell4.core.NumberObserver.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.NumberObserver.data" title="ecell4.core.NumberObserver.data"><code class="xref py py-obj docutils literal"><span class="pre">data</span></code></a></td>
<td>Return a list of the numbers of molecules you specified.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.NumberObserver.next_time" title="ecell4.core.NumberObserver.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the next time for logging.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.NumberObserver.num_steps" title="ecell4.core.NumberObserver.num_steps"><code class="xref py py-obj docutils literal"><span class="pre">num_steps</span></code></a></td>
<td>Return the number of steps.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.NumberObserver.reset" title="ecell4.core.NumberObserver.reset"><code class="xref py py-obj docutils literal"><span class="pre">reset</span></code></a></td>
<td>Reset the internal state.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.NumberObserver.targets" title="ecell4.core.NumberObserver.targets"><code class="xref py py-obj docutils literal"><span class="pre">targets</span></code></a></td>
<td>Return a list of <code class="docutils literal"><span class="pre">Species</span></code>, which this <code class="docutils literal"><span class="pre">Observer</span></code> observes</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.core.NumberObserver.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NumberObserver.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.NumberObserver.data">
<code class="descname">data</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.NumberObserver.data" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of the numbers of molecules you specified.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of lists of the number of molecules.</dt>
<dd>The size of a return value is equal to <code class="docutils literal"><span class="pre">num_steps</span></code>.
Each element of a return value is a list consisting of
time and the number of molecules specified at the construction.</dd>
</dl>
</dd>
</dl>
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
<dl class="docutils">
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of <code class="docutils literal"><span class="pre">Species</span></code>. This is generated from arguments</dt>
<dd>you gave at the construction.</dd>
</dl>
</dd>
</dl>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Observer">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Observer</code><a class="headerlink" href="#ecell4.core.Observer" title="Permalink to this definition">¶</a></dt>
<dd><p>A wrapper for a base class of Observers.</p>
<p>Warning: This is mainly for developers.
Do not use this for your simulation.</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Observer.next_time" title="ecell4.core.Observer.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the next time for logging.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Observer.reset" title="ecell4.core.Observer.reset"><code class="xref py py-obj docutils literal"><span class="pre">reset</span></code></a></td>
<td>Reset the internal state.</td>
</tr>
</tbody>
</table>
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
<dt id="ecell4.core.PlanarSurface">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">PlanarSurface</code><a class="headerlink" href="#ecell4.core.PlanarSurface" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a planar surface, which is available to define
structures.</p>
<p>PlanarSurface(origin, e0, e1)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.PlanarSurface.as_base" title="ecell4.core.PlanarSurface.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.PlanarSurface.dimension" title="ecell4.core.PlanarSurface.dimension"><code class="xref py py-obj docutils literal"><span class="pre">dimension</span></code></a></td>
<td>Return a dimension of this shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.PlanarSurface.is_inside" title="ecell4.core.PlanarSurface.is_inside"><code class="xref py py-obj docutils literal"><span class="pre">is_inside</span></code></a></td>
<td>Return if the given point is inside or not.</td>
</tr>
</tbody>
</table>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd>value (float): Zero or negative if the given point is inside.</dd>
</dl>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.ReactionRule">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">ReactionRule</code><a class="headerlink" href="#ecell4.core.ReactionRule" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a reaction rule between <code class="docutils literal"><span class="pre">Species</span></code>.</p>
<p>ReactionRule(reactants=None, products=None, k=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.ReactionRule.add_product" title="ecell4.core.ReactionRule.add_product"><code class="xref py py-obj docutils literal"><span class="pre">add_product</span></code></a></td>
<td>Append a product to the end.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.ReactionRule.add_reactant" title="ecell4.core.ReactionRule.add_reactant"><code class="xref py py-obj docutils literal"><span class="pre">add_reactant</span></code></a></td>
<td>Append a reactant to the end.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.ReactionRule.as_string" title="ecell4.core.ReactionRule.as_string"><code class="xref py py-obj docutils literal"><span class="pre">as_string</span></code></a></td>
<td>Return an unicode string describing this object.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.ReactionRule.count" title="ecell4.core.ReactionRule.count"><code class="xref py py-obj docutils literal"><span class="pre">count</span></code></a></td>
<td>Count the number of matches for reactants.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.ReactionRule.generate" title="ecell4.core.ReactionRule.generate"><code class="xref py py-obj docutils literal"><span class="pre">generate</span></code></a></td>
<td>Generate <a href="#id1"><span class="problematic" id="id2">``</span></a>ReactionRule``s from given reactants.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.ReactionRule.k" title="ecell4.core.ReactionRule.k"><code class="xref py py-obj docutils literal"><span class="pre">k</span></code></a></td>
<td>Return the kinetic rate constant as a float value.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.ReactionRule.products" title="ecell4.core.ReactionRule.products"><code class="xref py py-obj docutils literal"><span class="pre">products</span></code></a></td>
<td>List all products.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.ReactionRule.reactants" title="ecell4.core.ReactionRule.reactants"><code class="xref py py-obj docutils literal"><span class="pre">reactants</span></code></a></td>
<td>List all reactants.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.ReactionRule.set_k" title="ecell4.core.ReactionRule.set_k"><code class="xref py py-obj docutils literal"><span class="pre">set_k</span></code></a></td>
<td>Set a kinetic rate constant.</td>
</tr>
</tbody>
</table>
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
<dd><p>Generate <a href="#id3"><span class="problematic" id="id4">``</span></a>ReactionRule``s from given reactants.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd><dl class="first last docutils">
<dt>reactants (list): A list of <code class="docutils literal"><span class="pre">Species</span></code>. The order of <code class="docutils literal"><span class="pre">reactants</span></code></dt>
<dd>is respected.</dd>
</dl>
</dd>
<dt>Return:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of <a href="#id5"><span class="problematic" id="id6">``</span></a>ReactionRule``s. The reactants of each</dt>
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

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Rod">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Rod</code><a class="headerlink" href="#ecell4.core.Rod" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a Rod shape, which is available to define
structures. The cylinder is aligned to x-axis.</p>
<p>Rod(length, radius, origin=Real3(0, 0, 0))</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Rod.as_base" title="ecell4.core.Rod.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Rod.dimension" title="ecell4.core.Rod.dimension"><code class="xref py py-obj docutils literal"><span class="pre">dimension</span></code></a></td>
<td>Return a dimension of this shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Rod.distance" title="ecell4.core.Rod.distance"><code class="xref py py-obj docutils literal"><span class="pre">distance</span></code></a></td>
<td>Return a minimum distance from the given point to the surface.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Rod.is_inside" title="ecell4.core.Rod.is_inside"><code class="xref py py-obj docutils literal"><span class="pre">is_inside</span></code></a></td>
<td>Return if the given point is inside or not.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Rod.length" title="ecell4.core.Rod.length"><code class="xref py py-obj docutils literal"><span class="pre">length</span></code></a></td>
<td>Return a length of a cylinder part.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Rod.origin" title="ecell4.core.Rod.origin"><code class="xref py py-obj docutils literal"><span class="pre">origin</span></code></a></td>
<td>Return a center position of mass</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Rod.radius" title="ecell4.core.Rod.radius"><code class="xref py py-obj docutils literal"><span class="pre">radius</span></code></a></td>
<td>Return a radius of a cylinder.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Rod.shift" title="ecell4.core.Rod.shift"><code class="xref py py-obj docutils literal"><span class="pre">shift</span></code></a></td>
<td>Move the center toward the given displacement</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Rod.surface" title="ecell4.core.Rod.surface"><code class="xref py py-obj docutils literal"><span class="pre">surface</span></code></a></td>
<td>Create and return a surface shape.</td>
</tr>
</tbody>
</table>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>distance (float): A minimum distance from the given point.</dt>
<dd>Negative if the given point is inside.</dd>
</dl>
</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Rod.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd>value (float): Zero or negative if the given point is inside.</dd>
</dl>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>vec (Real3): A displacement.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Rod.surface">
<code class="descname">surface</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Rod.surface" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a surface shape.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd>shape (RodSurface): The surface shape.</dd>
</dl>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.RodSurface">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">RodSurface</code><a class="headerlink" href="#ecell4.core.RodSurface" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a hollow rod surface shape, which is
available to define structures. The cylinder is aligned to x-axis.</p>
<p>RodSurface(length, radius, origin=Real3(0, 0, 0))</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.RodSurface.as_base" title="ecell4.core.RodSurface.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.RodSurface.dimension" title="ecell4.core.RodSurface.dimension"><code class="xref py py-obj docutils literal"><span class="pre">dimension</span></code></a></td>
<td>Return a dimension of this shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.RodSurface.distance" title="ecell4.core.RodSurface.distance"><code class="xref py py-obj docutils literal"><span class="pre">distance</span></code></a></td>
<td>Return a minimum distance from the given point to the surface.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.RodSurface.inside" title="ecell4.core.RodSurface.inside"><code class="xref py py-obj docutils literal"><span class="pre">inside</span></code></a></td>
<td>Create and return a volume shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.RodSurface.is_inside" title="ecell4.core.RodSurface.is_inside"><code class="xref py py-obj docutils literal"><span class="pre">is_inside</span></code></a></td>
<td>Return if the given point is inside or not.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.RodSurface.length" title="ecell4.core.RodSurface.length"><code class="xref py py-obj docutils literal"><span class="pre">length</span></code></a></td>
<td>Return a length of a cylinder part.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.RodSurface.origin" title="ecell4.core.RodSurface.origin"><code class="xref py py-obj docutils literal"><span class="pre">origin</span></code></a></td>
<td>Return a center position of mass</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.RodSurface.radius" title="ecell4.core.RodSurface.radius"><code class="xref py py-obj docutils literal"><span class="pre">radius</span></code></a></td>
<td>Return a radius of a cylinder.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.RodSurface.shift" title="ecell4.core.RodSurface.shift"><code class="xref py py-obj docutils literal"><span class="pre">shift</span></code></a></td>
<td>Move the center toward the given displacement</td>
</tr>
</tbody>
</table>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>distance (float): A minimum distance from the given point.</dt>
<dd>Negative if the given point is inside.</dd>
</dl>
</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.RodSurface.inside">
<code class="descname">inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a volume shape.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd>shape (Rod): The volume shape.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.RodSurface.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.RodSurface.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd>value (float): Zero or negative if the given point is inside.</dd>
</dl>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>vec (Real3): A displacement.</dd>
</dl>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Shape">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Shape</code><a class="headerlink" href="#ecell4.core.Shape" title="Permalink to this definition">¶</a></dt>
<dd><p>A wrapper for a base class of Shapes.</p>
<p>Warning: This is mainly for developers.
Do not use this for your simulation.</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Shape.dimension" title="ecell4.core.Shape.dimension"><code class="xref py py-obj docutils literal"><span class="pre">dimension</span></code></a></td>
<td>Return a dimension of this shape.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Shape.is_inside" title="ecell4.core.Shape.is_inside"><code class="xref py py-obj docutils literal"><span class="pre">is_inside</span></code></a></td>
<td>Return if the given point is inside or not.</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.core.Shape.dimension">
<code class="descname">dimension</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Shape.dimension" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a dimension of this shape.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Shape.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Shape.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd>value (float): Zero or negative if the given point is inside.</dd>
</dl>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.Species">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Species</code><a class="headerlink" href="#ecell4.core.Species" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a type of molecules with attributes.</p>
<p>Species(serial=None, radius=None, D=None, location=None)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Species.add_unit" title="ecell4.core.Species.add_unit"><code class="xref py py-obj docutils literal"><span class="pre">add_unit</span></code></a></td>
<td>Append an <code class="docutils literal"><span class="pre">UnitSpecies</span></code> to the end.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Species.count" title="ecell4.core.Species.count"><code class="xref py py-obj docutils literal"><span class="pre">count</span></code></a></td>
<td>Count the number of matches for a pattern given as a <code class="docutils literal"><span class="pre">Species</span></code>.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Species.deserialize" title="ecell4.core.Species.deserialize"><code class="xref py py-obj docutils literal"><span class="pre">deserialize</span></code></a></td>
<td>Reset the serial.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Species.get_attribute" title="ecell4.core.Species.get_attribute"><code class="xref py py-obj docutils literal"><span class="pre">get_attribute</span></code></a></td>
<td>Return an attribute as an unicode string.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Species.has_attribute" title="ecell4.core.Species.has_attribute"><code class="xref py py-obj docutils literal"><span class="pre">has_attribute</span></code></a></td>
<td>Return if the attribute exists or not.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Species.list_attributes" title="ecell4.core.Species.list_attributes"><code class="xref py py-obj docutils literal"><span class="pre">list_attributes</span></code></a></td>
<td>List all attributes.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Species.num_units" title="ecell4.core.Species.num_units"><code class="xref py py-obj docutils literal"><span class="pre">num_units</span></code></a></td>
<td>Return the number of <code class="docutils literal"><span class="pre">UnitSpecies</span></code>.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Species.remove_attribute" title="ecell4.core.Species.remove_attribute"><code class="xref py py-obj docutils literal"><span class="pre">remove_attribute</span></code></a></td>
<td>Remove an attribute.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Species.serial" title="ecell4.core.Species.serial"><code class="xref py py-obj docutils literal"><span class="pre">serial</span></code></a></td>
<td>Return the serial name as an unicode string.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Species.set_attribute" title="ecell4.core.Species.set_attribute"><code class="xref py py-obj docutils literal"><span class="pre">set_attribute</span></code></a></td>
<td>Set an attribute.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Species.units" title="ecell4.core.Species.units"><code class="xref py py-obj docutils literal"><span class="pre">units</span></code></a></td>
<td>Return a list of all <code class="docutils literal"><span class="pre">UnitSpecies</span></code>.</td>
</tr>
</tbody>
</table>
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

<dl class="class">
<dt id="ecell4.core.Sphere">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">Sphere</code><a class="headerlink" href="#ecell4.core.Sphere" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a sphere shape, which is available to define
structures.</p>
<p>Sphere(center, radius)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Sphere.as_base" title="ecell4.core.Sphere.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Sphere.dimension" title="ecell4.core.Sphere.dimension"><code class="xref py py-obj docutils literal"><span class="pre">dimension</span></code></a></td>
<td>Return a dimension of this shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Sphere.distance" title="ecell4.core.Sphere.distance"><code class="xref py py-obj docutils literal"><span class="pre">distance</span></code></a></td>
<td>Return a minimum distance from the given point to the surface.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.Sphere.is_inside" title="ecell4.core.Sphere.is_inside"><code class="xref py py-obj docutils literal"><span class="pre">is_inside</span></code></a></td>
<td>Return if the given point is inside or not.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.Sphere.surface" title="ecell4.core.Sphere.surface"><code class="xref py py-obj docutils literal"><span class="pre">surface</span></code></a></td>
<td>Create and return a surface shape.</td>
</tr>
</tbody>
</table>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>distance (float): A minimum distance from the given point.</dt>
<dd>Negative if the given point is inside.</dd>
</dl>
</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Sphere.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Sphere.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd>value (float): Zero or negative if the given point is inside.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.Sphere.surface">
<code class="descname">surface</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.Sphere.surface" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a surface shape.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd>shape (SphericalSurface): The surface shape.</dd>
</dl>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.SphericalSurface">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">SphericalSurface</code><a class="headerlink" href="#ecell4.core.SphericalSurface" title="Permalink to this definition">¶</a></dt>
<dd><p>A class representing a hollow spherical surface, which is
available to define structures.</p>
<p>SphericalSurface(center, radius)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.SphericalSurface.as_base" title="ecell4.core.SphericalSurface.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.SphericalSurface.dimension" title="ecell4.core.SphericalSurface.dimension"><code class="xref py py-obj docutils literal"><span class="pre">dimension</span></code></a></td>
<td>Return a dimension of this shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.SphericalSurface.distance" title="ecell4.core.SphericalSurface.distance"><code class="xref py py-obj docutils literal"><span class="pre">distance</span></code></a></td>
<td>Return a minimum distance from the given point to the surface.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.SphericalSurface.inside" title="ecell4.core.SphericalSurface.inside"><code class="xref py py-obj docutils literal"><span class="pre">inside</span></code></a></td>
<td>Create and return a volume shape.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.SphericalSurface.is_inside" title="ecell4.core.SphericalSurface.is_inside"><code class="xref py py-obj docutils literal"><span class="pre">is_inside</span></code></a></td>
<td>Return if the given point is inside or not.</td>
</tr>
</tbody>
</table>
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
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>distance (float): A minimum distance from the given point.</dt>
<dd>Negative if the given point is inside.</dd>
</dl>
</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.SphericalSurface.inside">
<code class="descname">inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.SphericalSurface.inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Create and return a volume shape.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd>shape (Sphere): The volume shape.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.SphericalSurface.is_inside">
<code class="descname">is_inside</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.SphericalSurface.is_inside" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given point is inside or not.</p>
<dl class="docutils">
<dt>Args:</dt>
<dd>pos (Real3): A position.</dd>
<dt>Returns:</dt>
<dd>value (float): Zero or negative if the given point is inside.</dd>
</dl>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.core.TimeoutObserver">
<em class="property">class </em><code class="descclassname">ecell4.core.</code><code class="descname">TimeoutObserver</code><a class="headerlink" href="#ecell4.core.TimeoutObserver" title="Permalink to this definition">¶</a></dt>
<dd><p>An <a href="#id7"><span class="problematic" id="id8">``</span></a>Observer``class to stop simulation at the given calculation time.</p>
<p>TimeoutObserver(interval)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.TimeoutObserver.accumulation" title="ecell4.core.TimeoutObserver.accumulation"><code class="xref py py-obj docutils literal"><span class="pre">accumulation</span></code></a></td>
<td>Return the accumulation time.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.TimeoutObserver.as_base" title="ecell4.core.TimeoutObserver.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.TimeoutObserver.duration" title="ecell4.core.TimeoutObserver.duration"><code class="xref py py-obj docutils literal"><span class="pre">duration</span></code></a></td>
<td>Return the last time to be called.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.TimeoutObserver.interval" title="ecell4.core.TimeoutObserver.interval"><code class="xref py py-obj docutils literal"><span class="pre">interval</span></code></a></td>
<td>Return the timeout in seconds.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.TimeoutObserver.reset" title="ecell4.core.TimeoutObserver.reset"><code class="xref py py-obj docutils literal"><span class="pre">reset</span></code></a></td>
<td>Reset the internal state.</td>
</tr>
</tbody>
</table>
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
<dd><p>An <a href="#id9"><span class="problematic" id="id10">``</span></a>Observer``class to log the number of molecules just at the time
you assigned.</p>
<p>TimingNumberObserver(t, species)</p>
<p class="rubric">Methods</p>
<table border="1" class="longtable docutils">
<colgroup>
<col width="10%" />
<col width="90%" />
</colgroup>
<tbody valign="top">
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.TimingNumberObserver.as_base" title="ecell4.core.TimingNumberObserver.as_base"><code class="xref py py-obj docutils literal"><span class="pre">as_base</span></code></a></td>
<td>Clone self as a base class.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.TimingNumberObserver.data" title="ecell4.core.TimingNumberObserver.data"><code class="xref py py-obj docutils literal"><span class="pre">data</span></code></a></td>
<td>Return a list of the numbers of molecules you specified.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.TimingNumberObserver.next_time" title="ecell4.core.TimingNumberObserver.next_time"><code class="xref py py-obj docutils literal"><span class="pre">next_time</span></code></a></td>
<td>Return the next time for logging.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.TimingNumberObserver.num_steps" title="ecell4.core.TimingNumberObserver.num_steps"><code class="xref py py-obj docutils literal"><span class="pre">num_steps</span></code></a></td>
<td>Return the number of steps.</td>
</tr>
<tr class="row-odd"><td><a class="reference internal" href="#ecell4.core.TimingNumberObserver.reset" title="ecell4.core.TimingNumberObserver.reset"><code class="xref py py-obj docutils literal"><span class="pre">reset</span></code></a></td>
<td>Reset the internal state.</td>
</tr>
<tr class="row-even"><td><a class="reference internal" href="#ecell4.core.TimingNumberObserver.targets" title="ecell4.core.TimingNumberObserver.targets"><code class="xref py py-obj docutils literal"><span class="pre">targets</span></code></a></td>
<td>Return a list of <code class="docutils literal"><span class="pre">Species</span></code>, which this <code class="docutils literal"><span class="pre">Observer</span></code> observes</td>
</tr>
</tbody>
</table>
<dl class="method">
<dt id="ecell4.core.TimingNumberObserver.as_base">
<code class="descname">as_base</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimingNumberObserver.as_base" title="Permalink to this definition">¶</a></dt>
<dd><p>Clone self as a base class. This function is for developers.</p>
</dd></dl>

<dl class="method">
<dt id="ecell4.core.TimingNumberObserver.data">
<code class="descname">data</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.core.TimingNumberObserver.data" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of the numbers of molecules you specified.</p>
<dl class="docutils">
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of lists of the number of molecules.</dt>
<dd>The size of a return value is equal to <code class="docutils literal"><span class="pre">num_steps</span></code>.
Each element of a return value is a list consisting of
time and the number of molecules specified at the construction.</dd>
</dl>
</dd>
</dl>
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
<dl class="docutils">
<dt>Returns:</dt>
<dd><dl class="first last docutils">
<dt>list: A list of <code class="docutils literal"><span class="pre">Species</span></code>. This is generated from arguments</dt>
<dd>you gave at the construction.</dd>
</dl>
</dd>
</dl>
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
