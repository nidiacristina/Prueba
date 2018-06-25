#!/usr/bin/env python3

# takes an HDF5 file generated from the SMPM solver and produces an XDMF
# (v2) document representing its contents.
#
# the supplied HDF5 file's contents are inspected to determine what type
# of SMPM file it is (field output, initial conditions, restart file) and
# instantiate the class that represents its contents.  said instantiation
# knows how to generate tailored XDMF.

from __future__ import print_function

import getopt
import h5py
import numpy
import sys

def usage( script_name ):
    """
    Takes a name of the script (full path, name, etc) and prints its usage to standard
    output.
    """

    print( """Usage: {:s} [-h] <HDF5 file>

Takes an HDF5 file produced by the SMPM Incompressible Navier-Stokes solver and
generates a version 2 XDMF file on standard output that describes its contents.  The
XDMF file allows the HDF5 file to be visualized in VTK-based tools like ParaView,
Visit, and MayaVi.

The HDF5 file may be the output file containing field variables at each time
step, the initial conditions file used to start the solver, or the restart
file that is periodically written to checkpoint the solver's state.  The type
of file is automatically determined from the groups and datasets it contains
and unknown file types generate an error.

The command line options above are described below:

  -h        Displays this help message and exits.
""".format( script_name ) )

class SMPMFile( object ):
    """
    Class describing a generic SMPM-based HDF5 file.

    The XDMF generated contains a single computational domain with a
    structured, curvilinear mesh with a DOM like so:

      Xdmf
        Domain
          Topology
          Geometry

    Sub-classes of SMPMFile can reference the domain's grid either via
    XPath:

        <xi:include xpointer="element(/1/1/1)"/>
        <xi:include xpointer="element(/1/1/2)"/>

    Or via the Topology and Geometry nodes' Reference attribute:

        <Topology Reference="/Xdmf/Domain/Topology[1]"/>
        <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>

    """

    def __init__( self, h5, generators=None ):
        """
        Constructor that takes an HDF5 file and zero or more XDMF fragment
        generators.

        Takes 2 arguments:

          h5         - h5py.File object for the HDF5 file to generate XDMF
                       for.
          generators - Optional iterable that returns functions that produce
                       XDMF fragments as strings.  Each generator is invoked
                       sequentially and the output is appended to the XDMF
                       grid information.

        Returns 1 value:

          self - The newly constructed SMPMFile object.

        """

        # we start without any parameters.
        self._info = {}

        # XDMF generating functions used to fill content after the SMPM
        # solver's grid and before the XDMF footer.  iterable producing
        # functions that return XDMF fragments as strings.
        if generators is None:
            self._generators = []
        else:
            self._generators = generators

        # keep track of the file these parameters came from.
        self["file_name"] = h5.filename

        # store the basic grid parameters.
        self["n"]         = h5["/grid/n"].value
        self["mx"]        = h5["/grid/mx"].value * self["n"]
        self["mz"]        = h5["/grid/mz"].value * self["n"]

        # handle the y-dimension separately to be compatible with older
        # files and solver's that don't unconditionally write the third
        # dimension.
        #
        # XXX: remove this once the solver's restart capability supports
        #      3D.
        #
        if "my" in h5["/grid"]:
            self["my"] = h5["/grid/my"].value
        else:
            self["my"] = 1

        # work around applications that write scalars as 1x1 data sets.
        #
        # NOTE: this is needed to deal with Octave's 3rd party h5create()
        #       routine that, as of 2016/12/23, does not handle scalars.
        #       1x1 datasets are returned as NumPy ndarray's rather than
        #       a NumPy data type, so we take the first element of the
        #       array instead.
        #
        if isinstance( self["n"], numpy.ndarray ):
            self["n"] = self["n"][0]
        if isinstance( self["mx"], numpy.ndarray ):
            self["mx"] = self["mx"][0]
        if isinstance( self["my"], numpy.ndarray ):
            self["my"] = self["my"][0]
        if isinstance( self["mz"], numpy.ndarray ):
            self["mz"] = self["mz"][0]

    def __getitem__( self, key ):
        """
        Retrieves the value for an key within the parameters.

        Takes 1 argument:

          key - The key whose value is requested.

        Returns 1 value:

          value - The value associated with key.

        """

        return self._info[key]

    def __setitem__( self, key, value ):
        """
        Sets the value for an key within the parameters.

        Takes 2 arguments:

          key   - The key whose value will be set.
          value - The value to set.

        Returns nothing.

        """

        self._info[key] = value

    def xdmf_geometry_type( self ):
        """
        Gets the XDMF geometry type appropriate for the file.

        Takes no arguments.

        Returns 1 value:

          geometry - XDMF geometry string, suitable for use as the Geometry
                     node's GeometryType attribute's value.

        """

        #
        # NOTE: our grid coordinates are vectors of scalars rather than
        #       vectors of triplets, hence the use of X_Y(_Z) instead of
        #       XY(Z).
        #
        if self["my"] > 1:
            return "X_Y_Z"
        else:
            return "X_Y"

    def xdmf_grid_dimensions( self ):
        """
        Gets the dimensions of the XDMF grid in the order of what the SMPM
        solver writes out.

        Takes no arguments.

        Returns 1 value:

          grid_dimensions - XDMF dimensions, suitable for use as the DataItem
                            element's Dimension attribute's value.

        """

        # individual node locations are stored fastest to slowest (Z, X, then
        # Y AKA by planes), while the individual coordinates are described to match the
        # solver's view of the world (X-Z planes map to X-Y planes).
        grid_dimensions = "{mx:d} {mz:d}"

        #
        # NOTE: ParaView 5.1.2 (presumably other VTK-based visualization tools
        #       as well) does not handle a leading singleton dimension, so we
        #       only indicate data are 3D if there are three full dimensions.
        #
        if self["my"] > 1:
            grid_dimensions = "{my:d} " + grid_dimensions

        return grid_dimensions.format( mx=self["mx"],
                                       my=self["my"],
                                       mz=self["mz"] )

    def xdmf_topology_type( self ):
        """
        Gets the topoology type of the XDMF grid appropriate for the file.

        Takes no arguments.

        Returns 1 value:

          topology - XDMF topology, suitable for use as the Topology element's
                     TopologyType attribute's value.
        """

        if self["my"] > 1:
            return "3DSMesh"
        else:
            return "2DSMesh"

    def generate_xdmf( self ):
        """
        Constructs an XDMF document (v2) describing the file and returns it as a
        string.  The output of this object's generators, if any, is
        concatenated together and appended after the fragment describing the
        SMPM grid.

        Takes no arguments.

        Returns 1 value:

          xdmf - XDMF document (v2) serialized as a string.

        """

        xdmf = """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">
    <Domain>
{grid:s}
{content:s}
    </Domain>
</Xdmf>"""

        # call each of the generators in sequence and combine them into a
        # single, new-line delimited string.
        content = "\n".join( map( lambda x: x(), self._generators ) )

        return xdmf.format( grid=self.generate_xdmf_grid(),
                            content=content )

    def generate_xdmf_grid( self ):
        """
        Constructs an XDMF fragment describing the SMPM grid using a
        structured, curvilinear mesh where each coordinate is stored as a
        separate variable within the HDF5 described.

        Takes no arguments.

        Returns 1 value:

          grid - String containing the SMPM grid as an XDMF fragment.

        """

        # leave a comment for future archaeologists so they understand the
        # mapping between HDF5 data and the XDMF view of the world.
        #
        # NOTE: the extra indentation aligns the comment in the output.  this
        #       should be done in a better manner.
        grid_comment    = """Individual node locations are stored slowest to fastest (X then Z), while the
          order of individual coordinates match the solver's view of the world (the X-Z
          plane maps to the X-Y plane)."""
        geometry_type   = self.xdmf_geometry_type()
        grid_dimensions = self.xdmf_grid_dimensions()
        topology_type   = self.xdmf_topology_type()
        y_coord         = ""

        # a third dimension means we need to specify another coordinate as
        # well as explain why the SMPM solver's axes need to be permuted.
        if self["my"] > 1:
            grid_comment = """Individual node locations are stored slowest to fastest (Y, X, then Z,
          AKA by planes), while the order of individual coordinates match the solver's
          view of the world (X-Z planes map to X-Y planes)."""

            y_coord = """
            <DataItem Dimensions="{grid_dimensions:s}" NumberType="Float" Precision="8" Format="HDF">
            {file_name:s}:/grid/y
            </DataItem>""".format( file_name=self["file_name"],
                                   grid_dimensions=grid_dimensions )

        # NOTE: the geometry is specified X, Z, then Y to match XDMF's X, Y,
        #       and Z dimensions.
        grid = """        <!--
          {grid_comment:s}
        -->
        <Topology TopologyType="{topology_type:s}" NumberOfElements="{grid_dimensions:s}"/>
        <Geometry GeometryType="{geometry_type:s}">
            <DataItem Dimensions="{grid_dimensions:s}" NumberType="Float" Precision="8" Format="HDF">
            {file_name:s}:/grid/x
            </DataItem>
            <DataItem Dimensions="{grid_dimensions:s}" NumberType="Float" Precision="8" Format="HDF">
            {file_name:s}:/grid/z
            </DataItem>{y_coord:s}
        </Geometry>"""

        return grid.format( file_name=self["file_name"],
                            geometry_type=geometry_type,
                            grid_comment=grid_comment,
                            grid_dimensions=grid_dimensions,
                            topology_type=topology_type,
                            y_coord=y_coord )

class ICFile( SMPMFile ):
    """
    Class describing an initial conditions (IC) HDF5 file.
    """

    def __init__( self, h5 ):
        """
        Constructor that takes an initial conditions (IC) HDF5 file.

        Takes 1 argument:

          h5         - h5py.File object for the initial conditions file to
                       generate XDMF for.

        Returns 1 value:

          self - The newly constructed ICFile object.

        """

        super().__init__( h5, (self.generate_xdmf_initial_conditions,) )

    def generate_xdmf_initial_conditions( self ):
        """
        Constructs an XDMF fragment describing the SMPM initial conditions.
        Each of the field variables is described.

        Takes no arguments.

        Returns 1 value:

          ic - String containing the initial conditions as an XDMF fragment.

        """

        grid_dimensions = self.xdmf_grid_dimensions()
        uy              = ""

        if self["my"] > 1:
            uy = """
                <Attribute Name="uy" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/ic/uy
                    </DataItem>
                </Attribute>"""

        # XXX: uncomment the background current information when it is written
        #      as proper 3D data.
        ic = """
        <Grid Name="initial_conditions1" CollectionType="Spatial" GridType="Collection">
            <Grid Name="initial_conditions">

                <xi:include xpointer="element(/1/1/1)"/>
                <xi:include xpointer="element(/1/1/2)"/>

                <Attribute Name="ux" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/ic/ux
                    </DataItem>
                </Attribute>{uy:s}
                <Attribute Name="uz" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/ic/uz
                    </DataItem>
                </Attribute>
                <Attribute Name="rho" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/ic/rho
                    </DataItem>
                </Attribute>
                <Attribute Name="rho_bar" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/ic/rho_bar
                    </DataItem>
                </Attribute>
<!--
                <Attribute Name="ubc" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{mx:d} {mz:d}" Format="HDF">
                    {file_name:s}:/ic/ubc
                    </DataItem>
                </Attribute>
                <Attribute Name="dubcdz" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{mx:d} {mz:d}" Format="HDF">
                    {file_name:s}:/ic/dubcdz
                    </DataItem>
                </Attribute>
-->
            </Grid>
        </Grid>"""

        return ic.format( file_name=self["file_name"],
                          grid_dimensions=grid_dimensions,
                          # XXX: these are only needed for the 2D slice.
                          mx=self["mx"],
                          mz=self["mz"],
                          uy=uy.format( file_name=self["file_name"],
                                        grid_dimensions=grid_dimensions ) )

class OutputFile( SMPMFile ):
    """
    Class describing a field output HDF5 file.
    """
    def __init__( self, h5 ):
        """
        Constructor that takes a field output HDF5 file.

        Takes 1 argument:

          h5         - h5py.File object for the field output file to generate
                       XDMF for.

        Returns 1 value:

          self - The newly constructed OutputFile object.

        """

        super().__init__( h5, (self.generate_xdmf_timesteps,) )

        self["number_time_steps"] = h5['/field/number_steps'].value

    def generate_xdmf_timesteps( self ):
        """
        Constructs an XDMF fragment describing the SMPM solver's field outputs
        on a temporal grid.

        Takes no arguments.

        Returns 1 value:

          timesteps - String containing the field output as an XDMF fragment.

        """

        grid_dimensions = self.xdmf_grid_dimensions()
        uy              = ""
        timesteps       = """
        <Grid Name="mesh1" CollectionType="Temporal" GridType="Collection">
"""

        if self["my"] > 1:
            uy = """
                <Attribute Name="uy" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/field/step{time_step:d}/uy
                    </DataItem>
                </Attribute>"""

        # define a grid of scalar values (ux, uy (possibly), uz, and rho) for each
        # timestep written by the solver.  the topology and geometry points
        # back to the global versions.
        step_format_string = """
            <Grid Name="time{time_step:d}">
                <Time Value="{time_step:d}" />

                <xi:include xpointer="element(/1/1/1)"/>
                <xi:include xpointer="element(/1/1/2)"/>

                <Attribute Name="ux" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/field/step{time_step:d}/ux
                    </DataItem>
                </Attribute>{uy:s}
                <Attribute Name="uz" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/field/step{time_step:d}/uz
                    </DataItem>
                </Attribute>
                <Attribute Name="rho" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/field/step{time_step:d}/rho
                    </DataItem>
                </Attribute>{extra:s}
            </Grid>
"""

        ic_format_string = """
                <Attribute Name="rho_bar" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/field/step0/rho_bar
                    </DataItem>
                </Attribute>"""

        # write out the initial conditions as step #0.
        timesteps += step_format_string.format( time_step=0,
                                                extra=ic_format_string.format( file_name=self["file_name"],
                                                                               grid_dimensions=grid_dimensions ),
                                                file_name=self["file_name"],
                                                grid_dimensions=grid_dimensions,
                                                uy=uy.format( file_name=self["file_name"],
                                                              grid_dimensions=grid_dimensions,
                                                              time_step=0 ) )

        # write out the remaining steps.
        for step_number in range( self["number_time_steps"] ):
            timesteps += step_format_string.format( time_step=step_number + 1,
                                                    extra="",
                                                    file_name=self["file_name"],
                                                    grid_dimensions=grid_dimensions,
                                                    uy=uy.format( file_name=self["file_name"],
                                                                  grid_dimensions=grid_dimensions,
                                                                  time_step=step_number + 1 ) )

        timesteps += "        </Grid>"

        return timesteps

class RestartFile( SMPMFile ):
    """
    Class describing a restart HDF5 file.
    """
    def __init__( self, h5 ):
        """
        Constructor that takes a restart HDF5 file.

        Takes 1 argument:

          h5         - h5py.File object for the restart file to generate XDMF
                       for.

        Returns 1 value:

          self - The newly constructed RestartFile object.

        """
        super().__init__( h5, (self.generate_xdmf_restart_conditions,) )

    def generate_xdmf_restart_conditions( self ):
        """
        Constructs an XDMF fragment describing the SMPM restart checkpoint.
        Each of the field variables is described.

        Takes no arguments.

        Returns 1 value:

          restart - String containing the restart data as an XDMF fragment.

        """

        grid_dimensions = self.xdmf_grid_dimensions()

        uy  = ""
        Nuy = ""

        if self["my"] > 1:
            uy = """
                <Attribute Name="uy" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/uy
                    </DataItem>
                </Attribute>
                <Attribute Name="uy0" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/uy0
                    </DataItem>
                </Attribute>
                <Attribute Name="uy1" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/uy1
                    </DataItem>
                </Attribute>
                <Attribute Name="uy2" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/uy2
                    </DataItem>
                </Attribute>"""

            Nuy = """
                <Attribute Name="Nuy0" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nuy0
                    </DataItem>
                </Attribute>
                <Attribute Name="Nuy1" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nuy1
                    </DataItem>
                </Attribute>
                <Attribute Name="Nuy2" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nuy2
                    </DataItem>
                </Attribute>"""

        # XXX: uncomment the background current information when it is written
        #      as proper 3D data.
        restart = """
        <Grid Name="restart1" CollectionType="Spatial" GridType="Collection">
            <Grid Name="solver_conditions">

                <xi:include xpointer="element(/1/1/1)"/>
                <xi:include xpointer="element(/1/1/2)"/>

                <Attribute Name="ux" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/ux
                    </DataItem>
                </Attribute>
                <Attribute Name="ux0" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/ux0
                    </DataItem>
                </Attribute>
                <Attribute Name="ux1" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/ux1
                    </DataItem>
                </Attribute>
                <Attribute Name="ux2" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/ux2
                    </DataItem>
                </Attribute>{uy:s}
                <Attribute Name="uz" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/uz
                    </DataItem>
                </Attribute>
                <Attribute Name="uz0" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/uz0
                    </DataItem>
                </Attribute>
                <Attribute Name="uz1" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/uz1
                    </DataItem>
                </Attribute>
                <Attribute Name="uz2" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/uz2
                    </DataItem>
                </Attribute>
                <Attribute Name="rho" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/rho
                    </DataItem>
                </Attribute>
                <Attribute Name="rho0" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/rho0
                    </DataItem>
                </Attribute>
                <Attribute Name="rho1" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/rho1
                    </DataItem>
                </Attribute>
                <Attribute Name="rho2" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/rho2
                    </DataItem>
                </Attribute>
                <Attribute Name="rho_bar" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/rho_bar
                    </DataItem>
                </Attribute>
<!--
                <Attribute Name="ubc" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{mx:d} {mz:d}" Format="HDF">
                    {file_name:s}:/restart/ubc
                    </DataItem>
                </Attribute>
                <Attribute Name="dubcdz" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{mx:d} {mz:d}" Format="HDF">
                    {file_name:s}:/restart/dubcdz
                    </DataItem>
                </Attribute>
-->
                <Attribute Name="Nrho0" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nrho0
                    </DataItem>
                </Attribute>
                <Attribute Name="Nrho1" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nrho1
                    </DataItem>
                </Attribute>
                <Attribute Name="Nrho2" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nrho2
                    </DataItem>
                </Attribute>
                <Attribute Name="Nux0" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nux0
                    </DataItem>
                </Attribute>
                <Attribute Name="Nux1" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nux1
                    </DataItem>
                </Attribute>
                <Attribute Name="Nux2" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nux2
                    </DataItem>
                </Attribute>{Nuy:s}
                <Attribute Name="Nuz0" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nuz0
                    </DataItem>
                </Attribute>
                <Attribute Name="Nuz1" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nuz1
                    </DataItem>
                </Attribute>
                <Attribute Name="Nuz2" AttributeType="Scalar" Center="Node">
                    <DataItem DataType="Float" Precision="8" Dimensions="{grid_dimensions:s}" Format="HDF">
                    {file_name:s}:/restart/Nuz2
                    </DataItem>
                </Attribute>
            </Grid>
        </Grid>"""

        return restart.format( file_name=self["file_name"],
                               grid_dimensions=grid_dimensions,
                               # XXX: these are only needed for the 2D slice.
                               mx=self["mx"],
                               mz=self["mz"],
                               Nuy=Nuy,
                               uy=uy ).format( file_name=self["file_name"],
                                               grid_dimensions=grid_dimensions )

def main( argv ):
    """
    Takes a list of arguments, including the name of the application calling this method,
    and generates an XDMF file on standard output.
    """

    # parse our command line options.
    try:
        opts, args = getopt.getopt( argv[1:], "h" )
    except getopt.GetoptError as error:
        sys.stderr.write( "Error processing option: {:s}\n".format( str( error ) ) )
        sys.exit( 1 )

    # handle any valid options were were presented.
    for opt, arg in opts:
        if opt == '-h':
            usage( argv[0] )
            sys.exit()

    # make sure we were called properly.
    if len( args ) != 1:
        sys.stderr.write( "Expected 1 argument but received {:d}.\n".format( len( args ) ) )
        usage( argv[0] )
        sys.exit( 1 )

    # open the supplied file and pull out the variables we need to describe
    # its contents.
    h5 = h5py.File( args[0], "r" )

    # identify what kind of file we have from a unique dataset/group at the
    # root and construct the appropriate set of parameters.  unknown files
    # result in an error.
    mapping = { 'restart': RestartFile,
                'field':   OutputFile,
                'ic':      ICFile }
    for data_group in mapping:
        if data_group in h5:
            parms = mapping[data_group](h5)
            break
    else:
        sys.stderr.write( "'{:s}' does not appear to be a SMPM HDF5 file!\n".format( h5.filename ) )
        sys.exit( 1 )

    h5.close()

    # build XDMF based on the parameters supplied and dump it to standard
    # output.
    print( parms.generate_xdmf() )

if __name__ == "__main__":
    main( sys.argv )

