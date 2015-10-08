  <span class="target" id="module-ecell4.meso"></span><dl class="class">
<dt id="ecell4.meso.MesoscopicFactory">
<em class="property">class </em><code class="descclassname">ecell4.meso.</code><code class="descname">MesoscopicFactory</code><a class="headerlink" href="#ecell4.meso.MesoscopicFactory" title="Permalink to this definition">¶</a></dt>
<dd><p>A factory class creating a MesoscopicWorld instance and a MesoscopicSimulator instance.</p>
<p>MesoscopicFactory(matrix_sizes=None, GSLRandomNumberGenerator rng=None)</p>
<dl class="method">
<dt id="ecell4.meso.MesoscopicFactory.create_simulator">
<code class="descname">create_simulator</code><span class="sig-paren">(</span><em>arg1</em>, <em>arg2</em><span class="sig-paren">)</span> &rarr; MesoscopicSimulator<a class="headerlink" href="#ecell4.meso.MesoscopicFactory.create_simulator" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a MesoscopicSimulator instance.</p>
<dl class="docutils">
<dt>arg1 <span class="classifier-delimiter">:</span> <span class="classifier">MesoscopicWorld</span></dt>
<dd>a world</dd>
</dl>
<p>or</p>
<dl class="docutils">
<dt>arg1 <span class="classifier-delimiter">:</span> <span class="classifier">Model</span></dt>
<dd>a simulation model</dd>
<dt>arg2 <span class="classifier-delimiter">:</span> <span class="classifier">MesoscopicWorld</span></dt>
<dd>a world</dd>
</dl>
<dl class="docutils">
<dt>MesoscopicSimulator:</dt>
<dd>the created simulator</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicFactory.create_world">
<code class="descname">create_world</code><span class="sig-paren">(</span><em>arg1=None</em><span class="sig-paren">)</span> &rarr; MesoscopicWorld<a class="headerlink" href="#ecell4.meso.MesoscopicFactory.create_world" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a MesoscopicWorld instance.</p>
<dl class="docutils">
<dt>arg1 <span class="classifier-delimiter">:</span> <span class="classifier">Real3</span></dt>
<dd>The lengths of edges of a MesoscopicWorld created</dd>
</dl>
<p>or</p>
<dl class="docutils">
<dt>arg1 <span class="classifier-delimiter">:</span> <span class="classifier">str</span></dt>
<dd>The path of a HDF5 file for MesoscopicWorld</dd>
</dl>
<dl class="docutils">
<dt>MesoscopicWorld:</dt>
<dd>the created world</dd>
</dl>
</dd></dl>

</dd></dl>

<dl class="class">
<dt id="ecell4.meso.MesoscopicSimulator">
<em class="property">class </em><code class="descclassname">ecell4.meso.</code><code class="descname">MesoscopicSimulator</code><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator" title="Permalink to this definition">¶</a></dt>
<dd><p>A class running the simulation with the meso algorithm.</p>
<p>MesoscopicSimulator(m, w)</p>
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
<dl class="docutils">
<dt>list:</dt>
<dd>The list of reaction rules and infos.</dd>
</dl>
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
<dl class="docutils">
<dt>duration <span class="classifier-delimiter">:</span> <span class="classifier">Real</span></dt>
<dd>A duration for running a simulation.
A simulation is expected to be stopped at t() + duration.</dd>
<dt>observers <span class="classifier-delimiter">:</span> <span class="classifier">list of Obeservers, optional</span></dt>
<dd>Observers</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.set_dt">
<code class="descname">set_dt</code><span class="sig-paren">(</span><em>dt</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.set_dt" title="Permalink to this definition">¶</a></dt>
<dd><p>Set a step interval.</p>
<dl class="docutils">
<dt>dt <span class="classifier-delimiter">:</span> <span class="classifier">Real</span></dt>
<dd>A step interval</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.set_t" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the current time.</p>
<dl class="docutils">
<dt>t <span class="classifier-delimiter">:</span> <span class="classifier">Real</span></dt>
<dd>A current time.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicSimulator.step">
<code class="descname">step</code><span class="sig-paren">(</span><em>upto=None</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.meso.MesoscopicSimulator.step" title="Permalink to this definition">¶</a></dt>
<dd><p>Step the simulation.</p>
<dl class="docutils">
<dt>upto <span class="classifier-delimiter">:</span> <span class="classifier">Real, optional</span></dt>
<dd>The time which to step the simulation up to</dd>
</dl>
<dl class="docutils">
<dt>bool:</dt>
<dd>True if the simulation did not reach the given time.
When upto is not given, nothing will be returned.</dd>
</dl>
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
<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.add_molecules">
<code class="descname">add_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em>, <em>c=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.add_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Add some molecules.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>a species of molecules to add</dd>
<dt>num <span class="classifier-delimiter">:</span> <span class="classifier">Integer</span></dt>
<dd>the number of molecules to add</dd>
<dt>c <span class="classifier-delimiter">:</span> <span class="classifier">Integer or Shape, optional</span></dt>
<dd>a coordinate or shape to add molecules on</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.add_structure">
<code class="descname">add_structure</code><span class="sig-paren">(</span><em>sp</em>, <em>shape</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.add_structure" title="Permalink to this definition">¶</a></dt>
<dd><p>Add a structure.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>a species suggesting the shape.</dd>
<dt>shape <span class="classifier-delimiter">:</span> <span class="classifier">Shape</span></dt>
<dd>a shape of the structure.</dd>
</dl>
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
<dl class="docutils">
<dt>m <span class="classifier-delimiter">:</span> <span class="classifier">Model</span></dt>
<dd>a model to bind</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.check_structure">
<code class="descname">check_structure</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.check_structure" title="Permalink to this definition">¶</a></dt>
<dd><p>on_structure(sp, g) -&gt; bool</p>
<p>Check if the given subvolume is belonging to the structure.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>A species for the target structure.</dd>
<dt>g <span class="classifier-delimiter">:</span> <span class="classifier">Integer3</span></dt>
<dd>a global coordinate pointing a subvolume</dd>
</dl>
<dl class="docutils">
<dt>bool:</dt>
<dd>True if the subvolume is belonging to the structure.</dd>
</dl>
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
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>A species for the target structure.</dd>
<dt>g <span class="classifier-delimiter">:</span> <span class="classifier">Integer3</span></dt>
<dd>a global coordinate pointing a subvolume</dd>
</dl>
<dl class="docutils">
<dt>Real:</dt>
<dd>The occupancy of the structure.
As a default, return 1 if the subvolume overlaps with the structure,
and return 0 otherwise.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.get_volume">
<code class="descname">get_volume</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; Real<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.get_volume" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a volume of the given structure.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>A species for the target structure.</dd>
</dl>
<dl class="docutils">
<dt>Real:</dt>
<dd>A total volume of subvolumes belonging to the structure.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.global2coord">
<code class="descname">global2coord</code><span class="sig-paren">(</span><em>g</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.global2coord" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a global coordinate to a coordinate.</p>
<dl class="docutils">
<dt>g <span class="classifier-delimiter">:</span> <span class="classifier">Integer3</span></dt>
<dd>a global coordinate</dd>
</dl>
<dl class="docutils">
<dt>Integer:</dt>
<dd>A coordinate</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.has_structure">
<code class="descname">has_structure</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.has_structure" title="Permalink to this definition">¶</a></dt>
<dd><p>Return if the given structure is in the space or not.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>A species for the target structure.</dd>
</dl>
<dl class="docutils">
<dt>bool:</dt>
<dd>True if the given structure is in self.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.list_coordinates">
<code class="descname">list_coordinates</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [Integer]<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.list_coordinates" title="Permalink to this definition">¶</a></dt>
<dd><p>Return a list of coordinates of molecules belonging to the given species.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>A species of molecules.</dd>
</dl>
<dl class="docutils">
<dt>list:</dt>
<dd>A list of coordinates.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.list_coordinates_exact">
<code class="descname">list_coordinates_exact</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.list_coordinates_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>list_coordinates(sp) -&gt; [Integer]</p>
<p>Return a list of coordinates of molecules belonging to the given species.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>A species of molecules.</dd>
</dl>
<dl class="docutils">
<dt>list:</dt>
<dd>A list of coordinates.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.list_particles">
<code class="descname">list_particles</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.list_particles" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the list of particles.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species, optional</span></dt>
<dd>The species of particles to list up
If no species is given, return the whole list of particles.</dd>
</dl>
<dl class="docutils">
<dt>list:</dt>
<dd>The list of particles (of the given species)</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.list_particles_exact">
<code class="descname">list_particles_exact</code><span class="sig-paren">(</span><em>sp</em><span class="sig-paren">)</span> &rarr; [(ParticleID, Particle)]<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.list_particles_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the list of particles of a given species.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>The species of particles to list up</dd>
</dl>
<dl class="docutils">
<dt>list:</dt>
<dd>The list of particles of a given species</dd>
</dl>
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
<dl class="docutils">
<dt>filename <span class="classifier-delimiter">:</span> <span class="classifier">str</span></dt>
<dd>a filename to load from</dd>
</dl>
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
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>A species whose molecules you count</dd>
<dt>c <span class="classifier-delimiter">:</span> <span class="classifier">Integer, optional</span></dt>
<dd>A coordinate.</dd>
</dl>
<dl class="docutils">
<dt>Integer:</dt>
<dd>The number of molecules (of a given species)</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.num_molecules_exact">
<code class="descname">num_molecules_exact</code><span class="sig-paren">(</span><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.num_molecules_exact" title="Permalink to this definition">¶</a></dt>
<dd><p>num_particles_exact(sp, c=None) -&gt; Integer</p>
<p>Return the number of molecules within the suggested subvolume.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>The species of molecules to count</dd>
<dt>c <span class="classifier-delimiter">:</span> <span class="classifier">Integer, optional</span></dt>
<dd>A coordinate.</dd>
</dl>
<dl class="docutils">
<dt>Integer:</dt>
<dd>The number of molecules of a given species</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.num_subvolumes">
<code class="descname">num_subvolumes</code><span class="sig-paren">(</span><em>sp=None</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.num_subvolumes" title="Permalink to this definition">¶</a></dt>
<dd><p>Return the number of subvolumes.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species, optional</span></dt>
<dd>A species specifying a structure.
When no species is given, return the total number of subvolumes.</dd>
</dl>
<dl class="docutils">
<dt>Integer:</dt>
<dd>The number of subvolumes.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.on_structure">
<code class="descname">on_structure</code><span class="sig-paren">(</span><em>sp</em>, <em>g</em><span class="sig-paren">)</span> &rarr; bool<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.on_structure" title="Permalink to this definition">¶</a></dt>
<dd><p>Check if the given species would be on the proper structure at the coordinate.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>a species scheduled to be placed</dd>
<dt>g <span class="classifier-delimiter">:</span> <span class="classifier">Integer3</span></dt>
<dd>a global coordinate pointing a subvolume</dd>
</dl>
<dl class="docutils">
<dt>bool:</dt>
<dd>True if it is on the proper structure.</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.position2coordinate">
<code class="descname">position2coordinate</code><span class="sig-paren">(</span><em>pos</em><span class="sig-paren">)</span> &rarr; Integer<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.position2coordinate" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a position to a coordinate.</p>
<dl class="docutils">
<dt>pos <span class="classifier-delimiter">:</span> <span class="classifier">Real3</span></dt>
<dd>A position</dd>
</dl>
<dl class="docutils">
<dt>Integer:</dt>
<dd>A coordinate</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.position2global">
<code class="descname">position2global</code><span class="sig-paren">(</span><em>pos</em><span class="sig-paren">)</span> &rarr; Integer3<a class="headerlink" href="#ecell4.meso.MesoscopicWorld.position2global" title="Permalink to this definition">¶</a></dt>
<dd><p>Transform a position to a global coordinate.</p>
<dl class="docutils">
<dt>pos <span class="classifier-delimiter">:</span> <span class="classifier">Real3</span></dt>
<dd>A position</dd>
</dl>
<dl class="docutils">
<dt>Integer3:</dt>
<dd>A global coordinate</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.remove_molecules">
<code class="descname">remove_molecules</code><span class="sig-paren">(</span><em>sp</em>, <em>num</em>, <em>c=None</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.remove_molecules" title="Permalink to this definition">¶</a></dt>
<dd><p>Remove the molecules.</p>
<dl class="docutils">
<dt>sp <span class="classifier-delimiter">:</span> <span class="classifier">Species</span></dt>
<dd>a species whose molecules to remove</dd>
<dt>num <span class="classifier-delimiter">:</span> <span class="classifier">Integer</span></dt>
<dd>a number of molecules to be removed</dd>
<dt>c <span class="classifier-delimiter">:</span> <span class="classifier">Integer, optional</span></dt>
<dd>A coordinate.</dd>
</dl>
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
<dl class="docutils">
<dt>filename <span class="classifier-delimiter">:</span> <span class="classifier">str</span></dt>
<dd>a filename to save to</dd>
</dl>
</dd></dl>

<dl class="method">
<dt id="ecell4.meso.MesoscopicWorld.set_t">
<code class="descname">set_t</code><span class="sig-paren">(</span><em>t</em><span class="sig-paren">)</span><a class="headerlink" href="#ecell4.meso.MesoscopicWorld.set_t" title="Permalink to this definition">¶</a></dt>
<dd><p>Set the value of the time of the world.</p>
<dl class="docutils">
<dt>t <span class="classifier-delimiter">:</span> <span class="classifier">Real</span></dt>
<dd>The time of the world</dd>
</dl>
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
