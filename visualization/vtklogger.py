import os
import shutil
import numpy
from vtk_xml_serial_unstructured import *
from _gfrd import Sphere, Cylinder, Box, Plane
from multi import Multi
from utils import crossproduct, INF


class VTKLogger:
    """
    """
    def __init__(self, sim, dir='vtkdata', buffer_size=None, show_shells=True, 
                 extra_particle_step=True, color_dict=None):
        """Create a new VTKLogger to visualize a simulation with 
        VTK/ParaView.

        Arguments:
            - sim
                an EGFRDSimulator.
            - dir
                the directory to which the data will be written. Can be 
                nested as in 'data/vtkdata'. By default the directory 
                'vtkdata' wil be used.
            - buffer_size
                an integer representing the maximum number of 
                simulation steps that will be stored in a 
                first-in-first-out buffer. The content of the buffer 
                will not be written to files until you call the 'stop' 
                method. By default no buffer is used, and simulation 
                data is written to file immediately after you call the 
                'log' method.
            - show_shells
                a boolean value which indicates whether besides 
                particle data also shell data should be recorded. True 
                by default.
            - extra_particle_step
                a boolean value which indicates whether besides each 
                simulation step an extra step should be recorded where 
                only the active particle is updated (it's shell stays 
                unchanged). True by default, which makes an eGFRD 
                simulations easier to understand.
            - color_dict
                a Python dictionary with species serial as a key and 
                color index as a value. The color index will later be 
                used to lookup the color for the particle of this 
                species in a lookuptable. You can choose which 
                lookuptable to use from within ParaView.
                The first Species you added to the model will have 
                serial 1, the second Species serial 2, etc.
                If you do not specify a color dictionary, the first 
                Species will also get color index 1, the second Species 
                color index 2, etc.

        To visualize the cylinders a workaround using tensors is used, as 
        explained here:
        http://www.paraview.org/pipermail/paraview/2009-March/011256.html.
        The mentioned tensorGlyph.xml should be supplied with this package.

        """
        self.sim = sim
        self.show_shells = show_shells
        self.extra_particle_step = extra_particle_step
        self.color_dict = color_dict

        self.vtk_writer = VTK_XML_Serial_Unstructured()

        # Filename stuff.
        assert len(dir) > 0
        self.dir = dir

        files_directory = self.dir + '/files'
        if os.path.exists(files_directory):
            shutil.rmtree(files_directory)
        os.makedirs(files_directory) 

        self.file_list = []
        # Static list of files used for regions and surfaces.
        self.static_list = []

        # First in first out buffer.
        self.buffer_size = buffer_size
        self.buffer = []

        # Step counter.
        self.i = 0

        # Needed for ParaView time hack. Note: don't make delta_t too 
        # small, should be relative to max time.
        self.delta_t = 1e-11
        self.last_time = INF

    def log(self):
        """Write the particle and shell data (last one is optional) for 
        the current timestep to .vtk files or to a buffer.

        """
        time = self.sim.t
        if (abs(time - self.last_time) < 1e-9) and self.show_shells:
            # Hack to make ParaView understand this is a different 
            # simulator step but with the same time.
            #
            # This is needed because:
            # 1. During multi steps the global time is not updated.
            # 2. During initialization the time is 0.
            # 3. Sometimes shells have mobility_radius = 0 --> dt = 0.
            #
            # Every step should be recorded, even when nothing has 
            # changed, so it is easier to relate an event index from 
            # the terminal output of a simulation with a visualization 
            # step in ParaView (divide by 2 actually).
            #
            # Now ParaView should perceive a state change.
            # Only needed when showing shells (this is a choice).
            time = self.last_time + self.delta_t

        # Get data.
        particle_data = self.get_particle_data()
        sphere_data, cylinder_data = self.get_sphere_and_cylinder_data()

        # Write to buffer or file.
        if self.buffer_size:
            # Store in buffer, instead of writing to file directly.
            self.buffer.append((time, self.i, particle_data, sphere_data, 
                                cylinder_data))

            if self.i > self.buffer_size:
                # FIFO.
                # Store 1 step more than defined buffer size: situation 
                # before, n steps, situation after.
                del self.buffer[0]
        else:
            # Write normal log.
            self.writelog(time, self.i,
                          (particle_data, sphere_data, cylinder_data))

        self.i += 1
        self.last_time = time

    def writelog(self, time, index,
                 (particle_data, sphere_data, cylinder_data)):
        # There are 4 cases to consider (a bit tricky):
        # I.   i=0, extra_particle_step=True
        # II.  i=1, extra_particle_step=True.
        # III. i=0, extra_particle_step=False.
        # IV.  i=1, extra_particle_step=False.

        if self.show_shells and self.extra_particle_step:
            # Show step where only particle_data have been updated.
            if index == 0:
                # Case I.
                sphere_data_1 = self.get_dummy_sphere_data()
                cylinder_data_1 = self.get_dummy_cylinder_data()
            else:
                # Case II.
                sphere_data_1 = self.previous_sphere_data
                cylinder_data_1 = self.previous_cylinder_data

            index *= 2;

            self.make_snapshot('particle_data', particle_data, index, time)
            self.make_snapshot('sphere_data', sphere_data_1, index, time)
            self.make_snapshot('cylinder_data', cylinder_data_1, index, time)

            # Continue with case. 
            time += self.delta_t / 2
            index += 1

        if index == 0:
            # Case III.
            sphere_data_2 = self.get_dummy_sphere_data()
            cylinder_data_2 = self.get_dummy_cylinder_data()
        else:
            # Case IV.
            sphere_data_2 = sphere_data
            cylinder_data_2 = cylinder_data

        self.make_snapshot('particle_data', particle_data, index, time)
        self.make_snapshot('sphere_data', sphere_data_2, index, time)
        self.make_snapshot('cylinder_data', cylinder_data_2, index, time)

        self.previous_sphere_data = sphere_data
        self.previous_cylinder_data = cylinder_data

    def make_snapshot(self, type, data, index='', time=None):
        # Write data to .vtu files.

        doc = self.vtk_writer.create_doc(data)
        file_name = 'files/' + type + str(index) + '.vtu'

        self.vtk_writer.write_doc(doc, self.dir + '/' + file_name)

        # Store filename and time in file_list, used by 
        # vtk_writer.write_pvd().
        if time == None:
            # This is a region or surface.
            self.static_list.append((type, file_name, None, None))
        else:
            self.file_list.append((type, file_name, index, time))

    def stop(self):
        """Write a list with all the .vtu files to a .pvd file.

        """
        # In case user forgot this, doesn't hurt.
        self.log()

        # Write contents of buffer.
        for index, entry in enumerate(self.buffer):
            # Use enumerate because ParaView wants index to start at 0; 
            # we can not use realindex.
            # An entry contains: (time, realindex, particle_data, 
            # sphere_data, cylinder_data)
            if index == 0:
                # Add dummy data with color range to 
                # particle/sphere/cylinder data.
                #
                # Particle/sphere/cylinder data contain 4 lists:
                #   - positions
                #   - radii
                #   - colors
                #   - tensors
                #
                dummy_particle_data = self.get_dummy_particle_data()
                particle_data = entry[2]
                particle_data[0].extend(dummy_particle_data[0])
                particle_data[1].extend(dummy_particle_data[1])
                particle_data[2].extend(dummy_particle_data[2])
                particle_data[3].extend(dummy_particle_data[3])

                dummy_sphere_data = self.get_dummy_sphere_data()
                sphere_data = entry[3]
                sphere_data[0].extend(dummy_sphere_data[0])
                sphere_data[1].extend(dummy_sphere_data[1])
                sphere_data[2].extend(dummy_sphere_data[2])
                sphere_data[3].extend(dummy_sphere_data[3])

                dummy_cylinder_data = self.get_dummy_cylinder_data()
                cylinder_data = entry[4]
                cylinder_data[0].extend(dummy_cylinder_data[0])
                cylinder_data[1].extend(dummy_cylinder_data[1])
                cylinder_data[2].extend(dummy_cylinder_data[2])
                cylinder_data[3].extend(dummy_cylinder_data[3])

            if index % 100 == 0 and index > 0:
                print('vtklogger finished writing from buffer up to step %s' % 
                      index)

            self.writelog(entry[0], index, entry[2:])

        # Write data for regions and surfaces only once.
        self.make_snapshot('cylindrical_surfaces', 
                           self.get_cylindrical_surface_data())
        self.make_snapshot('planar_surfaces',
                           self.get_planar_surface_data())
        self.make_snapshot('cuboidal_regions', 
                           self.get_cuboidal_region_data())

        # Finally, write PVD files.
        self.vtk_writer.write_pvd(self.dir + '/' + 'files.pvd', 
                                  self.file_list)

        self.vtk_writer.write_pvd(self.dir + '/' + 'static.pvd', 
                                  self.static_list)

    def get_particle_data(self):
        particles, colors = [], []

        # Get particle data from simulator.
        for species in self.sim.world.species:
            for particle_id in self.sim.world.get_particle_ids(species):
                particle = self.sim.world.get_particle(particle_id)[1]
                particles.append(particle)
                try:
                    color = self.color_dict[species.id.serial]
                except:
                    color = species.id.serial
                colors.append(color)

        if self.i == 0:
            dummy_particles, dummy_colors = \
                self.get_dummy_particles_and_color_range()
            particles.extend(dummy_particles)
            colors.extend(dummy_colors)

        return self.process_spheres(particles, colors)

    def get_sphere_and_cylinder_data(self):
        spheres, sphere_colors = [], []
        cylinders, cylinder_colors = [], []

        if self.show_shells == True:
            number_of_shells = self.sim.scheduler.size
        else:
            number_of_shells = 0

        if number_of_shells > 0:
            top_event_id = self.sim.scheduler.top[0]

        for object in self.sim.domains.itervalues():
            # Single and pairs.
            color = object.multiplicity

            # Highlight top_event for singles and pairs.
            if object.event_id == top_event_id:
                color = 0

            # Multi.
            if isinstance(object, Multi):
                color = 3

            try:
                shell = object.shell_list[0][1]
                # Only cylinders have size.
                shell.shape.size
                cylinders.append(shell)
                cylinder_colors.append(color)
            except:
                # Spheres: single, pair or multi.
                for _, shell in object.shell_list:
                    spheres.append(shell)
                    sphere_colors.append(color)

        return self.process_spheres(spheres, sphere_colors), \
               self.process_cylinders(cylinders, cylinder_colors)

    def get_cuboidal_region_data(self):
        boxes = [self.sim.world.get_structure("world").shape]

        return self.process_boxes(boxes)

    def get_planar_surface_data(self):
        world = self.sim.world.get_structure("world")
        boxes = [surface.shape for surface
                               in self.sim.world.structures
                               if isinstance(surface.shape, Plane)
                               and not surface == world]
    
        return self.process_boxes(boxes)

    def get_cylindrical_surface_data(self):
        # Todo. Make DNA blink when reaction takes place.
        cylinders = [surface for surface
                             in self.sim.world.structures
                             if isinstance(surface.shape, Cylinder)]
        return self.process_cylinders(cylinders)

    def process_spheres(self, spheres=[], color_list=[]):
        # Return 4 lists:
        #   - positions
        #   - radii
        #   - colors
        #   - tensors (not used)
        #
        pos_list, radius_list = [], []
        if len(spheres) == 0:
            # Add dummy sphere to stop ParaView from complaining.
            spheres = [self.get_dummy_sphere()]
            color_list = [0]

        for sphere in spheres:
            # Todo. Which one is it.
            try:
                position = sphere.position
                radius = sphere.radius
            except:
                position = sphere.shape.position
                radius = sphere.shape.radius
            self.append_lists(pos_list, position, radius_list, radius)

        return (pos_list, radius_list, color_list, [])

    def process_cylinders(self, cylinders=[], color_list=[]):
        # Return 4 lists:
        #   - positions
        #   - radii (not used)
        #   - colors
        #   - tensors
        #
        pos_list, tensor_list = [], []

        if len(cylinders) == 0:
            # Add dummy cylinder to stop TensorGlyph from complaining.
            cylinders = [self.get_dummy_cylinder()]
            color_list = [0]

        for cylinder in cylinders:
            # Todo. Which one is it.
            try:
                position = cylinder.position
                radius = cylinder.radius
                orientation = cylinder.unit_z
                size = cylinder.size
            except:
                position = cylinder.shape.position
                radius = cylinder.shape.radius
                orientation = cylinder.shape.unit_z
                size = cylinder.shape.size

            # Construct tensor. Use TensorGlyph plugin from:
            # http://www.paraview.org/pipermail/paraview/2009-March/011256.html
            # Unset Extract eigenvalues.

            # Select basis vector in which orientation is smallest.
            _, basis_vector = min(zip(abs(orientation), [[1, 0, 0], 
                                                         [0, 1, 0],
                                                         [0, 0, 1]]))
            # Find 2 vectors perpendicular to orientation.
            perpendicular1 = numpy.cross(orientation, basis_vector)
            perpendicular2 = numpy.cross(orientation, perpendicular1)
            # A 'tensor' is represented as an array of 9 values.
            # Stupid ParaView wants  a normal vector to the cylinder to 
            # orient it. So orientation and perpendicular1 swapped.
            tensor = numpy.concatenate((perpendicular1 * radius, 
                                        orientation * size,
                                        perpendicular2 * radius))

            self.append_lists(pos_list, position, tensor_list=tensor_list, 
                              tensor=tensor)

        return (pos_list, [], color_list, tensor_list)

    def process_boxes(self, boxes=[], color_list=[]):
        pos_list, tensor_list = [], []

        if len(boxes) == 0:
            # Add dummy box to stop TensorGlyph from complaining.
            boxes = [self.get_dummy_box()]
            color_list = [0]

        for box in boxes:
            try:
                dz = box.unit_z * box.extent[2]
            except AttributeError:
                # Planes don't have z dimension.
                unit_z = crossproduct(box.unit_x, box.unit_y)
                dz = unit_z * 1e-20
            tensor = numpy.concatenate((box.unit_x * box.extent[0],
                                        box.unit_y * box.extent[1],
                                        dz))
            self.append_lists(pos_list, box.position, tensor_list=tensor_list, 
                              tensor=tensor)

        return (pos_list, [], color_list, tensor_list)

    def append_lists(self, pos_list, pos, radius_list=[], radius=None, 
                     tensor_list=[], tensor=None):
        # Helper method.

        factor = 1
        pos_list.append(pos * factor)

        # Multiply radii and tensors by 2 because ParaView sets radius 
        # to 0.5 by default, and it wants full lengths for cylinders 
        # and we are storing half lengths.
        if radius:
            radius_list.append(radius * 2 * factor)

        if tensor != None:
            tensor = tensor * 2 * factor
            tensor_list.append(tensor)
    
    def get_dummy_particle_data(self):
        dummy_particles, dummy_colors = \
            self.get_dummy_particles_and_color_range()
        return self.process_spheres(dummy_particles, dummy_colors)

    def get_dummy_particles_and_color_range(self):
        # Start with dummy particle with color 1
        # (species.id.serial starts at 1).
        particles = [self.get_dummy_sphere()]
        colors = [1]

        # Add dummy particle with color is highest species index.
        particles.append(self.get_dummy_sphere())
        try:
            max_color = max(self.color_dict.values())
        except:
            for species in self.sim.world.species:
                # Find species with highest id.
                pass
            max_color = species.id.serial
        colors.append(max_color)

        return particles, colors

    def get_dummy_sphere_data(self):
        # Start with dummy sphere with color is 0.
        spheres = [self.get_dummy_sphere()]
        colors = [0]

        # Add dummy sphere with color is 3.
        spheres.append(self.get_dummy_sphere())
        colors.append(3)

        return self.process_spheres(spheres, colors)

    def get_dummy_cylinder_data(self):
        # Start with dummy cylinder with color is 0.
        cylinders = [self.get_dummy_cylinder()]
        colors = [0]

        # Add dummy cylinder with color is 3.
        cylinders.append(self.get_dummy_cylinder())
        colors.append(3)

        return self.process_cylinders(cylinders, colors)

    def get_dummy_sphere(self):
        return Sphere([0, 0, 0], 1e-20)

    def get_dummy_cylinder(self):
        return Cylinder([0, 0, 0], 1e-20, [0, 0, 1], 1e-20)

    def get_dummy_box(self):
        return Box([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
                   1e-20, 1e-20, 1e-20)

