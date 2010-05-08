#!/usr/bin/env python
import xml.dom.minidom
#import xml.dom.ext # python 2.5 and later


# Adapted from: 
# http://www.shocksolution.com/microfluidics-and-biotechnology/python-vtk-paraview/


class VTK_XML_Serial_Unstructured:
    def __init__(self):
        pass

    def create_doc(self, (pos_list, radii, colors, tensors)):
        # Document and root element
        doc = xml.dom.minidom.Document()
        root_element = doc.createElementNS("VTK", "VTKFile")
        root_element.setAttribute("type", "UnstructuredGrid")
        root_element.setAttribute("version", "0.1")
        root_element.setAttribute("byte_order", "LittleEndian")
        doc.appendChild(root_element)

        # Unstructured grid element
        unstructured_grid = doc.createElementNS("VTK", "UnstructuredGrid")
        root_element.appendChild(unstructured_grid)

        # Piece 0 (only one)
        # The "Piece" elements are meant for multiple pieces of 
        # *geometry*.  They are meant for streaming computation to 
        # reduce memory usage.  All the pieces have to have the same 
        # set of data arrays.
        # So we can not use that to group particles, spheres and 
        # cylinders into 1 file. Use .pvd file using parts for that.
        piece = doc.createElementNS("VTK", "Piece")
        piece.setAttribute("NumberOfPoints", str(len(pos_list)))
        piece.setAttribute("NumberOfCells", "0")
        unstructured_grid.appendChild(piece)


        # Points.
        points = doc.createElementNS("VTK", "Points")
        piece.appendChild(points)

        # Point location data.
        point_coords = doc.createElementNS("VTK", "DataArray")
        point_coords.setAttribute("type", "Float32")
        point_coords.setAttribute("format", "ascii")
        point_coords.setAttribute("NumberOfComponents", "3")
        points.appendChild(point_coords)

        string = str()
        for pos in pos_list:
            string = string + str(pos[0]) + ' ' + str(pos[1]) + \
                     ' ' + str(pos[2]) + ' '
        point_coords_data = doc.createTextNode(string)
        point_coords.appendChild(point_coords_data)


        # Cells.
        # Don't remove.
        cells = doc.createElementNS("VTK", "Cells")
        piece.appendChild(cells)

        # Cell locations.
        cell_connectivity = doc.createElementNS("VTK", "DataArray")
        cell_connectivity.setAttribute("type", "Int32")
        cell_connectivity.setAttribute("Name", "connectivity")
        cell_connectivity.setAttribute("format", "ascii")        
        cells.appendChild(cell_connectivity)

        # Cell location data.
        connectivity = doc.createTextNode("0")
        cell_connectivity.appendChild(connectivity)

        cell_offsets = doc.createElementNS("VTK", "DataArray")
        cell_offsets.setAttribute("type", "Int32")
        cell_offsets.setAttribute("Name", "offsets")
        cell_offsets.setAttribute("format", "ascii")                
        cells.appendChild(cell_offsets)
        offsets = doc.createTextNode("0")
        cell_offsets.appendChild(offsets)

        cell_types = doc.createElementNS("VTK", "DataArray")
        cell_types.setAttribute("type", "UInt8")
        cell_types.setAttribute("Name", "types")
        cell_types.setAttribute("format", "ascii")                
        cells.appendChild(cell_types)
        types = doc.createTextNode("1")
        cell_types.appendChild(types)


        # Point data.
        point_data = doc.createElementNS("VTK", "PointData")
        piece.appendChild(point_data)

        # Radii.
        if len(radii) > 0:
            radii_node = doc.createElementNS("VTK", "DataArray")
            radii_node.setAttribute("Name", "radii")
            radii_node.setAttribute("type", "Float32")
            radii_node.setAttribute("format", "ascii")
            point_data.appendChild(radii_node)

            string = str()
            for radius in radii:
                string = string + str(radius) + ' '
            radii_data = doc.createTextNode(string)
            radii_node.appendChild(radii_data)

        # Colors.
        if len(colors) > 0:
            color_node = doc.createElementNS("VTK", "DataArray")
            color_node.setAttribute("Name", "colors")
            color_node.setAttribute("NumberOfComponents", "1")
            color_node.setAttribute("type", "Float32")
            color_node.setAttribute("format", "ascii")
            point_data.appendChild(color_node)

            string = str()
            for color in colors:
                string = string + str(color) + ' '
            color_data = doc.createTextNode(string)
            color_node.appendChild(color_data)

        # Tensors.
        if len(tensors) > 0:

            # Hack to make VTK understand I want to color the 
            # TensorGlyphs.

            # I think there is a bug actually in VTK somewhere, when 
            # using tensorGlyph.xml to make a vtkTensorGlyph object 
            # with both a Tensor array and a Scalar array. When you 
            # select a Tensor or a Scalar array from the dropdown menu 
            # and click 'apply', SetInputArrayToProcess is called with 
            # idx = 0 both times. This is wrong.
            # 1. First element Tensor array gets overwritten.
            # 2. Scalar value is never written (which is accessed using 
            # GetInputArrayToProcess with an idx of 1).
            
            # The workaround here uses an additional Vector array (
            # vtk_tensor_glyph doesn't have a vector array), which when 
            # it gets updated actually results in 
            # SetInputArrayToProcess to be called with idx = 1. So by 
            # also supplying Paraview with a vector for each color 
            # (just 3 times the color int), the Tensor array doesn't 
            # get overwritten and GetInputArrayToProcess with idx of 1 
            # is also happy and retrieves the correct scalars.

            # TensorGlyph:
            # Tensors: 0
            # Scalars: 1 (but writes to 0 due to bug)
            # Vectors: 2 (but writes to 1 due to bug)

            # Glyph:
            # Scalars: 0
            # Vectors: 1

            # Tried to find the root cause using GDB but failed.

            # From a Python script you can use 
            # SelectInputScalars/Vectors/Tensors with the appropriate 
            # idx as first argument and all is fine. Not specifying a 
            # number ('') is equal to '0'.

            # Important files:
            #    ParaView3/VTK/Graphics/vtkTensorGlyph.cxx
            #    ParaView3/VTK/Filtering/vtkAlgorithm.cxx

            if len(colors) > 0:
                color_node = doc.createElementNS("VTK", "DataArray")
                color_node.setAttribute("Name", "colors_as_vectors")
                color_node.setAttribute("NumberOfComponents", "3")
                color_node.setAttribute("type", "Float32")
                color_node.setAttribute("format", "ascii")
                point_data.appendChild(color_node)

                string = str()
                for color in colors:
                    string += str(color) + ' ' + str(color) + ' ' + \
                              str(color) + ' '
                color_data = doc.createTextNode(string)
                color_node.appendChild(color_data)


            tensor_node = doc.createElementNS("VTK", "DataArray")
            tensor_node.setAttribute("Name", "tensors")
            tensor_node.setAttribute("NumberOfComponents", "9")
            tensor_node.setAttribute("type", "Float32")
            tensor_node.setAttribute("format", "ascii")
            point_data.appendChild(tensor_node)

            string = str()
            for tensor in tensors:
                for value in tensor:
                    # A 'tensor' is represented as a list of 9 values.
                    string = string + str(value) + ' '
            tensor_data = doc.createTextNode(string)
            tensor_node.appendChild(tensor_data)

        # Cell data (Dummy).
        cell_data = doc.createElementNS("VTK", "CellData")
        piece.appendChild(cell_data)

        return doc

    def write_doc(self, doc, file_name):
        # Write to file and exit.
        out_file = open(file_name, 'w')
        # xml.dom.ext.PrettyPrint(doc, file)
        doc.writexml(out_file, newl='\n')
        out_file.close()

    def write_pvd(self, file, file_list):
        out_file = open(file, 'w')

        pvd = xml.dom.minidom.Document()
        pvd_root = pvd.createElementNS("VTK", "VTKFile")
        pvd_root.setAttribute("type", "Collection")
        pvd_root.setAttribute("version", "0.1")
        pvd_root.setAttribute("byte_order", "LittleEndian")
        pvd.appendChild(pvd_root)

        collection = pvd.createElementNS("VTK", "Collection")
        pvd_root.appendChild(collection)

        for type, file_name, index, time in file_list:
            data_set = pvd.createElementNS("VTK", "DataSet")
            if time != None:
                # Problem with adding real time can be that 
                # TimestepValues is not updated in Proxy group="misc" 
                # in .pvsm file after a reload.
                time = str(time)
                data_set.setAttribute("timestep", time)
            data_set.setAttribute("group", "")
            data_set.setAttribute("part", type) # Use Extractblock.
            data_set.setAttribute("file", file_name)
            collection.appendChild(data_set)

        out_file = open(file, 'w')
        pvd.writexml(out_file, newl='\n')
        out_file.close()

