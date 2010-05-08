from paraview import servermanager
import os

READERS   = None
PARTICLES = None
SPHERES   = None
CYLINDERS = None
HELIX     = None
SURFACES  = None
DARK_BACKGROUND = None

# Settings
READERS   = True
PARTICLES = True
SPHERES   = True
CYLINDERS = True
HELIX     = True
SURFACES  = True

DARK_BACKGROUND = True

PARTICLE_RADIUS_SCALE_FACTOR = 1
HELIX_RADIUS_SCALE_FACTOR = 1
RESOLUTION = 18


def load_xml_plugin(filename):
    # Adapted from LoadPlugin in simple.py (ParaView 3.6).
    f = open(filename, 'r')
    if os.name == "posix":
        import libvtkPVServerManagerPython as libvtk
        parser = libvtk.vtkSMXMLParser()
    else:
        import vtkPVServerCommonPython as libvtk
        parser = libvtk.vtkSMXMLParser()
    if not parser.Parse(f.read()):
        raise RuntimeError("Problem loading plugin %s: %s" %
                           filename)
    parser.ProcessConfiguration(
            libvtk.vtkSMObject.GetProxyManager())
    # Update the modules
    servermanager.updateModules()

def MakeLT(rgbPoints):
    # Helper function.

    lt = servermanager.rendering.PVLookupTable()
    servermanager.Register(lt)
    lt.RGBPoints = rgbPoints

    return lt

def MakeBlueToRedLT(min, max):
    # Adapted from MakeBlueToRedLT in simple.py (ParaView 3.6).
    rgbPoints = [min, 0, 0, 1, max, 1, 0, 0]

    lt = MakeLT(rgbPoints)
    lt.ColorSpace = "HSV"

    return lt

def MakeCoolToWarmLT(min, max):
    r=0.137255
    g=0.239216
    b=0.709804
    r2=0.67451
    g2=0.141176
    b2=0.12549
    rgbPoints = [min, r, g, b, max, r2, g2, b2]

    lt = MakeLT(rgbPoints)
    lt.ColorSpace = "Diverging"

    return lt

def MakeNiceLT(min, max):
    # Define RGB points. These are tuples of 4 values. First one is
    # the scalar values, the other 3 the RGB values. 

    top_color = [1, 0, 0] # red
    single_color = [80. / 256, 160. / 256, 1] # light blue
    pair_color = [40. / 256, 80. / 256, 1] # dark blue
    multi_color = [160. / 256, 80. / 256, 1]

    rgbPoints = []
    rgbPoints.append(0)
    rgbPoints.extend(top_color)

    rgbPoints.append(1)
    rgbPoints.extend(single_color)

    rgbPoints.append(2)
    rgbPoints.extend(pair_color)

    rgbPoints.append(3)
    rgbPoints.extend(multi_color)
    
    return MakeLT(rgbPoints)

class _funcs_internals:
    # Helper class. Taken from simple.py (ParaView 3.6).
    "Internal class."
    first_render = True
    view_counter = 0
    rep_counter = 0

class Pipeline(object):
    """
    """
    def __init__(self, egfrd_directory, simulation_data_directory):
        """Create, initialize and build a new ParaView pipeline that 
        visualizes eGFRD simulation data produced by VTKLogger.

        Call this from the Python shell inside ParaView. Tested with 
        ParaView 3.4, 3.6 and 3.8. You might see some warnings, but it 
        should generally work.

        Arguments:
            - egfrd_directory
                path of the egfrd directory, for example: /home/user/egfrd
            - simulation_data_directory
                path of the directory that contains the files 
                'files.pvd' and 'static.pvd', as created by VTKLogger.

        """
        self.paraview_scripts_directory = egfrd_directory + '/visualization/'
        self.simulation_data_directory = simulation_data_directory
        
        self.initialize()
        self.rebuild()

    def initialize(self):
        dir = self.paraview_scripts_directory

        tensor_glyph_path = dir + '/tensorGlyph.xml'
        tensor_glyph_custom_path = dir + '/tensorGlyphWithCustomSource.xml'
        assert os.path.isfile(tensor_glyph_path), \
                'File does not exist: ' + tensor_glyph_path

        assert os.path.isfile(tensor_glyph_custom_path), \
                'File does not exist: ' + tensor_glyph_custom_path

        if not servermanager.ActiveConnection:
            exit('pvpython not supported. Use ParaView\'s Python shell.')

        # Detect ParaView version and load tensor plugins.
        try:
            from paraview import simple
            # This is ParaView 3.6 or higher.
            version = 6
            simple.LoadPlugin(tensor_glyph_path)
            simple.LoadPlugin(tensor_glyph_custom_path)
        except ImportError:
            # This must be ParaView 3.4.
            simple = None
            version = 4

            load_xml_plugin(tensor_glyph_path)
            load_xml_plugin(tensor_glyph_custom_path)

        # The RenderView is needed by clear() as well as by build().
        if len(servermanager.GetRenderViews()) > 0:
            rv = servermanager.GetRenderViews()[0]
        else:
            rv = servermanager.CreateRenderView()

        # Disable polygon offset (z-fighting prevention), otherwise 
        # wireframe is shifted.
        # http://www.mail-archive.com/paraview@paraview.org/msg08812.html
        from libvtkRenderingPython import vtkMapper
        vtkMapper.SetResolveCoincidentTopologyToOff()

        self.version = version
        self.rv = rv
        self.simple = simple

    def clear(self):
        version = self.version
        rv = self.rv
        simple = self.simple

        print 'Clear the pipeline.'

        # Reset time so that color range is detected correctly on build().
        rv.ViewTime = 0
        rv.StillRender()

        def name(proxy):
            # Return name of proxy.
            return (type(proxy)).__name__

        def cmp_tubes_filters_glyphs_blocks(x,y):
            # Using this function to sort the proxies will assure they are 
            # removed in the right order.
            if name(x) in ['GenerateTubes', 'TubeFilter', 'Tube']:
                return -1
            elif name(y) in ['GenerateTubes', 'TubeFilter', 'Tube']:
                return 1
            if name(x) == 'ProgrammableFilter':
                return -1
            elif name(y) == 'ProgrammableFilter':
                return 1
            elif name(x) == 'Glyph' or name(x)[:11] == 'TensorGlyph':
                return -1
            elif name(y) == 'Glyph' or name(y)[:11] == 'TensorGlyph':
                return 1
            if name(x) == 'ExtractBlock':
                return -1
            elif name(y) == 'ExtractBlock':
                return 1
            return cmp(x,y)

        # Remove lookup tables first.
        pxm = servermanager.ProxyManager()
        for proxy in pxm.GetProxiesInGroup('lookup_tables').itervalues():
            servermanager.UnRegister(proxy)

        if version == 4:
            # Then remove the source proxies.
            for proxy in sorted(pxm.GetProxiesInGroup('sources').itervalues(),
                                cmp_tubes_filters_glyphs_blocks):
                if name(proxy) == 'TensorGlyphWithCustomSource':
                    # Do nothing.
                    # Setting Source or Input gives: 
                    # 'QAbstractItemModel::endRemoveRows:
                    # Invalid index ( 2 , 0 ) in model 
                    # pqPipelineModel(0x26340b0)'
                    # http://www.paraview.org/Bug/view.php?id=9312
                    pass
                else:
                    # Avoid error:
                    # 'Connection sink not found in the pipeline model'.
                    if hasattr(proxy, "Source"):
                        proxy.Source = None
                    if hasattr(proxy, "Input"):
                        proxy.Input = None

                servermanager.UnRegister(proxy)

            # Finally remove the representations.
            for proxy in pxm.GetProxiesInGroup('representations').itervalues():
                servermanager.UnRegister(proxy)

            rv.Representations = []

        else:
            for proxy in sorted(simple.GetSources().itervalues(), 
                                cmp_tubes_filters_glyphs_blocks):
                # Avoid error:
                # 'Connection sink not found in the pipeline model'.
                if hasattr(proxy, "Input"):
                    proxy.Input = None
                if hasattr(proxy, "GlyphType"):
                    proxy.GlyphType = None

                simple.Delete(proxy)

        rv.ResetCamera()
        rv.StillRender()

    def rebuild(self):
        """First clear the pipeline, and then rebuild it.

        Use this method when your simulation data has changed.

        """
        version = self.version
        rv = self.rv
        data_dir = self.simulation_data_directory
        scripts_dir = self.paraview_scripts_directory

        self.clear()
        print 'Build the pipeline.'

        if READERS:
            files_pvd_path = data_dir + '/files.pvd'
            static_pvd_path = data_dir + '/static.pvd'
            assert os.path.isfile(files_pvd_path), \
                   'File does not exist: ' + files_pvd_path
            assert os.path.isfile(static_pvd_path), \
                   'File does not exist: ' + static_pvd_path

            files = self.add_pvd_reader(files_pvd_path, 'files.pvd')

            static = self.add_pvd_reader(static_pvd_path, 'static.pvd')


        if PARTICLES:
            particle_data = self.add_extract_block(files, [2], 'b1')

            particles = self.add_sphere_glyph(particle_data, name='Particles')
            particles.SetScaleFactor = PARTICLE_RADIUS_SCALE_FACTOR

            rep1 = self.show(particles)
            self.set_color(particles, rep1)


        if SPHERES:
            sphere_data = self.add_extract_block(files, [4], 'b2')

            spheres = self.add_sphere_glyph(sphere_data, RESOLUTION, 
                                            name='Spheres')

            rep2 = self.show(spheres)
            self.set_color(spheres, rep2, color_map=MakeNiceLT)
            rep2.Representation = 'Wireframe'
            rep2.Opacity = 0.25


        if CYLINDERS:
            cylinder_data = self.add_extract_block(files, [6], 'b3')

            cylinders = self.add_tensor_glyph(cylinder_data, 'Cylinder', 
                                              resolution=RESOLUTION, name='tg')

            programmable_filter = self.add_color_hack(cylinders,
                                                      name='Cylinders')
            rep3 = self.show(programmable_filter)
            self.set_color(programmable_filter, rep3, color_map=MakeNiceLT)
            rep3.Representation = 'Wireframe'
            rep3.Opacity = 1.0


        if SURFACES:
            # Cylindrical surfaces.
            cylindrical_surface_data = self.add_extract_block(static, [2], 'b4')


            if not HELIX:
                cylindrical_surfaces = \
                    self.add_tensor_glyph(cylindrical_surface_data, 'Cylinder',
                                          name='Cylindrical Surfaces',
                                          scale=HELIX_RADIUS_SCALE_FACTOR) 

                rep4 = self.show(cylindrical_surfaces)
                rep4.Representation = 'Wireframe'
                rep4.Opacity = 0.5
            else:
                helix_path = scripts_dir + '/helix.py'
                assert os.path.isfile(helix_path), \
                        'File does not exist: ' + helix_path
                helix_file = open(helix_path, 'r')
                helix_source = servermanager.sources.ProgrammableSource()
                helix_source.Script = 'HELIX_RADIUS_SCALE_FACTOR = ' + \
                                      str(HELIX_RADIUS_SCALE_FACTOR) + \
                                      '\n' + helix_file.read()
                helix_source.UpdatePipeline()
                servermanager.Register(helix_source, registrationName='ps')

                helix_file.close()


                tensor_glyph = self.add_tensor_glyph_with_custom_source(
                        cylindrical_surface_data, helix_source, name='tgwcs')
                if version == 4:
                    double_helix = \
                        servermanager.filters.TubeFilter(Input=tensor_glyph)
                else:
                    try:
                        double_helix = servermanager.filters.GenerateTubes(     
                                Input=tensor_glyph)
                    except AttributeError:
                        # ParaView 3.8.
                        double_helix = \
                            servermanager.filters.Tube(Input=tensor_glyph)
                servermanager.Register(double_helix,
                                       registrationName='Double Helix')

                # Compute helix radius.
                di = cylindrical_surface_data.GetDataInformation()
                pdi = di.GetPointDataInformation()
                tensor = pdi.GetArrayInformation('tensors')

                cylindrical_surface_radius = 1e100
                for i in range(9):
                    # Find minimum value of tensor larger than 0.
                    value = tensor.GetComponentRange(i)[0]
                    if value > 0 and value < cylindrical_surface_radius:
                        cylindrical_surface_radius = value

                helix_radius = HELIX_RADIUS_SCALE_FACTOR * \
                               cylindrical_surface_radius

                # Make double_helix a bit thinner than helix.
                double_helix.Radius = helix_radius / 20

                rep4 = self.show(double_helix)
                self.set_color(double_helix, rep4,
                               color_array_name = 'TubeNormals')


            # Planar surfaces.
            planar_surface_data = self.add_extract_block(static, [4], 'b5')
            planar_surfaces = self.add_tensor_glyph(planar_surface_data, 'Box', 
                                                    name='Planar Surfaces')

            rep5 = self.show(planar_surfaces)
            rep5.Representation = 'Surface'
            rep5.Opacity = 0.5
            rep5.Opacity = 0.5


            # Cuboidal surfaces.
            cuboidal_region_data = self.add_extract_block(static, [6], 'b6')
            cuboidal_regions = self.add_tensor_glyph(cuboidal_region_data,
                                                     'Box',
                                                     name='Cuboidal Regions')

            rep6 = self.show(cuboidal_regions)
            rep6.Representation = 'Wireframe'
            rep6.Opacity = 1.0
            rep6.LineWidth = 2 #1

            if DARK_BACKGROUND:
                color = [1, 1, 1] # White.
            else:
                color = [0, 0, 0] # Black.

            if version == 4:
                rep6.Color = color
            else:
                rep6.AmbientColor = color


        # Set camera.
        cam = rv.GetActiveCamera()
        # Sets focalpoint to center of box.
        rv.ResetCamera()
        cam.SetViewUp(0,0,1)
        focal = cam.GetFocalPoint()
        # Straigh in front of box.
        cam.SetPosition(focal[0]*10, focal[1], focal[2])

        if DARK_BACKGROUND:
            rv.Background = [0, 0, 0] # Black.
        else:
            rv.Background = [1, 1, 1] # White.

        # Changes cam.GetPosition() only by zooming in/out.
        rv.ResetCamera()
        rv.StillRender()

    def add_pvd_reader(self, file, name):
        # Helper function.
        version = self.version

        reader = servermanager.sources.PVDReader(FileName=file)
        servermanager.Register(reader, registrationName=name)
        if version == 4:
            reader.UpdatePipeline()
            # Update TimestepValues.
            reader.UpdatePipelineInformation();
        else:
            pass

        return reader

    def add_extract_block(self, data, indices, name):
        # Helper function.

        block = servermanager.filters.ExtractBlock(Input=data, 
                                                   BlockIndices=indices)

        # This is needed to make SetScaleFactor and TensorGlyph work.
        block.UpdatePipeline();

        servermanager.Register(block, registrationName=name)
        return block

    def add_sphere_glyph(self, input, resolution=None, name=None):
        # Helper function.
        version = self.version

        if version == 4:
            source = servermanager.sources.SphereSource()
            if resolution != None:
                source.ThetaResolution = resolution
                source.PhiResolution = resolution
            glyph = servermanager.filters.Glyph(Input=input, 
                                                Source=source)

            # Prevent "selected proxy value not in the list".
            # public.kitware.com/pipermail/paraview/2008-March/007416.html
            sourceProperty = glyph.GetProperty("Source")
            domain = sourceProperty.GetDomain("proxy_list");
            domain.AddProxy(source.SMProxy)

            glyph.SetScaleMode = 0
            # Look at the documentation of 
            # vtkAlgorithm::SetInputArrayToProcess() for details.
            glyph.SelectInputScalars = ['0', '', '', '', 'radii']

        else:
            glyph = servermanager.filters.Glyph(Input=input, 
                                                GlyphType='Sphere')
            glyph.ScaleMode = 'scalar'

            if resolution != None:
                glyph.GlyphType.ThetaResolution = resolution
                glyph.GlyphType.PhiResolution = resolution

        glyph.SetScaleFactor = 1

        if name != None:
            servermanager.Register(glyph, registrationName=name)
        else:
            servermanager.Register(glyph)

        return glyph

    def add_tensor_glyph(self, input, type, resolution=None, name=None, 
                         scale=None):
        # Helper function.
        version = self.version

        if version == 4:
            if type == 'Cylinder':
                source = servermanager.sources.CylinderSource()

                if resolution != None:
                    source.Resolution = resolution

                if scale != None:
                    source.Radius = 0.5 * scale

            elif type == 'Box':
                source = servermanager.sources.CubeSource()

            tensor_glyph = servermanager.filters.TensorGlyph(Input=input,
                                                             Source=source)

            tensor_glyph.SelectInputTensors = ['0', '', '', '', 'tensors']

            # The specified scalar array is the only array that gets copied.
            tensor_glyph.SelectInputScalars = ['1', '', '', '', 'colors']
        else:
            tensor_glyph = servermanager.filters.TensorGlyph(Input=input,
                                                             GlyphType=type)

            # Heads up. The first or the specified vector array is the only 
            # array that gets copied (scalar arrays don't get copied).
            tensor_glyph.Vectors = ['POINTS', 'colors_as_vectors']

            if resolution != None:
                tensor_glyph.GlyphType.Resolution = resolution

            if scale != None:
                tensor_glyph.GlyphType.Radius *= scale

        if name != None:
            servermanager.Register(tensor_glyph, registrationName=name)
        else:
            servermanager.Register(tensor_glyph)

        return tensor_glyph

    def add_tensor_glyph_with_custom_source(self, input, source, name=None):
        # Helper function.
        version = self.version

        Type = servermanager.filters.TensorGlyphWithCustomSource
        if version == 4:
            tensor_glyph = Type(Input=input, Source=source)
            tensor_glyph.SelectInputTensors = ['0', '', '', '', 'tensors']
        else:
            tensor_glyph = Type(Input=input, GlyphType=source)

        if name != None:
            servermanager.Register(tensor_glyph, registrationName=name)
        else:
            servermanager.Register(tensor_glyph)

        return tensor_glyph

    def add_color_hack(self, tensor_glyph, name):
        # Helper function.

        # Dealing with composite datasets:
        # http://www.itk.org/Wiki/Python_Programmable_Filter
        #
        # http://www.paraview.org/pipermail/paraview/2009-March/011267.html
        # As explainded in the above link, to give cylinders the right 
        # color, there are 2 options.
        # 1. Build ParaView from source, after adding in file
        # Paraview3/VTK/Graphics/vtkTensorGlyph.cxx after
        #   newScalars = vtkFloatArray::New();
        # a new line:
        #   newScalars->SetName("colors");
        #
        # 2. Use this Python script as a programmable filter.
        filter = servermanager.filters.ProgrammableFilter()
        filter.Initialize()
        filter.Input = tensor_glyph
        filter.Script = """
def flatten(input, output):
    output.ShallowCopy(input) # DeepCopy doesn't copy to tensor_glyph.
    output.GetPointData().GetArray(0).SetName('colors')
    # GetScalars() doesn't work in ParaView 3.4.
    #output.GetPointData().GetScalars().SetName('colors')

input = self.GetInput()
output = self.GetOutput()

output.CopyStructure(input)
iter = input.NewIterator()
iter.UnRegister(None)
iter.InitTraversal()
while not iter.IsDoneWithTraversal():
    curInput = iter.GetCurrentDataObject()
    curOutput = curInput.NewInstance()
    curOutput.UnRegister(None)
    output.SetDataSet(iter, curOutput)
    flatten(curInput, curOutput)
    iter.GoToNextItem()
"""

        filter.UpdatePipeline()

        servermanager.Register(filter, registrationName=name)

        return filter

    def set_color(self, proxy, rep, color_array_name='colors', 
                  color_map=MakeBlueToRedLT):
        # Helper function.
        version = self.version
        simple = self.simple

        rep.ColorArrayName = color_array_name

        if version == 4:
            rep.ColorAttributeType = 0 # point data

            pdi = proxy.GetDataInformation().GetPointDataInformation()
            color_array = pdi.GetArrayInformation(color_array_name)
            range = color_array.GetComponentRange(0)
        else:
            range = proxy.PointData.GetArray(color_array_name).GetRange()

        rep.LookupTable = color_map(range[0], range[1])

    def show(self, proxy):
        # Helper function.
        version = self.version
        rv = self.rv
        simple = self.simple

        if version == 4:
            # Adapted from Show in simple.py (ParaView 3.6).
            rep = servermanager.CreateRepresentation(proxy, rv)
            servermanager.ProxyManager().RegisterProxy("representations",
              "my_representation%d" % _funcs_internals.rep_counter, rep)
            _funcs_internals.rep_counter += 1
        else:
            rep = simple.Show(proxy)

        rep.Visibility = 1

        return rep

