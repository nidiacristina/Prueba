subroutine write_xdmf_field_file( xdmf_file_name, hdf5_file_name )
! Writes an XDMF file that represents the contents of the time-evolving HDF5
! field file written by the solver.  The file generated, along with the field
! file, may be visualized by anything that understands XDMF data including
! Paraview, VisIt, Mayavi, anything VTK-based, etc.  Creates field variable
! references for a fixed number of steps which includes one additional step to
! represent the field's initial conditions.
!
! NOTE: The XML written is not validated before it is generated and may
!       be invalid if care is not taken.
!
! NOTE: We jump through hoops to describe the output data as 2D when the
!       number of transverse planes is 1.  This is done to avoid problems with
!       VTK-based visualization tools (specifically, Paraview 5.1.2) which
!       crash when a 3D volume's leading dimension is 1.  Should this not be a
!       VTK issue and rather Greg's inability to understand the intracies of
!       XDMF v2, VTK, and ParaView then the mess below can be greatly
!       simplified.
!

  implicit none

  character(len=*), intent(in)        :: xdmf_file_name
  character(len=*), intent(in)        :: hdf5_file_name

  integer, parameter                  :: unit = 1024

  open( unit, file=xdmf_file_name )

  call write_xdmf_header( unit )

  ! Write out the grid and its topology.
  call write_xdmf_file_geometry( unit, hdf5_file_name )

  write( unit, '(A)' ) '        <Grid Name="mesh1" CollectionType="Temporal" GridType="Collection">'

  call write_xdmf_output_timestep( unit, hdf5_file_name )

  call write_xdmf_footer( unit )

  close( unit )

end subroutine write_xdmf_field_file

subroutine write_xdmf_output_timestep( unit, hdf5_file_name )
! Writes out field variables for a single timestep in XDMF.  This correctly
! handles writing out the initial conditions as step #0.

  use constants, only: n, nsubz, nsuby, nsubx
  use io_field, only:  field_name_initial_conditions_rho_bar, field_name_initial_conditions_ubc, &
                       field_name_initial_conditions_rho_bar_z, field_name_initial_conditions_dubcdz, &
                       field_name_rho, &
                       field_name_field, field_name_ux, field_name_uy, field_name_uz, &
                       field_name_x, field_name_y, field_name_z

  implicit none

  integer, intent(in)          :: unit
  character(len=*), intent(in) :: hdf5_file_name

  character(len=32)            :: dimensions_string

  call get_xdmf_dimensions_string( nsubx * n, nsuby, nsubz * n, dimensions_string )

  ! Write out a single timestep that inherits the topology of the grid.  Each
  ! field variable (rho, ux, uy, and uz) is referenced for all timesteps,
  ! while the background density stratification (rho_bar) is referenced for
  ! the initial conditions timestep (step #0).

  write( unit, '(A)' )              ''
  write( unit, '(A,I0,A)' )         '            <Grid Name="field">'
  write( unit, '(A)' )              ''
  write( unit, '(A)' )              '                <xi:include xpointer="element(/1/1/1)"/>'
  write( unit, '(A)' )              '                <xi:include xpointer="element(/1/1/2)"/>'
  write( unit, '(A)' )              ''
  write( unit, '(A)' )              '                <Attribute Name="ux" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A,A,A)' )    '                    ', hdf5_file_name, ':',  field_name_field, '/', field_name_ux
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  if ( nsuby > 1 ) then
     write( unit, '(A)' )           '                <Attribute Name="uy" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )       '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A,A,A)' ) '                    ', hdf5_file_name, ':', field_name_field, '/', field_name_uy
     write( unit, '(A)' )           '                    </DataItem>'
     write( unit, '(A)' )              '                </Attribute>'
  end if
  write( unit, '(A)' )              '                <Attribute Name="uz" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A,A,A)' )    '                    ', hdf5_file_name, ':', field_name_field, '/', field_name_uz
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="rho" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A,A,A)' )    '                    ', hdf5_file_name, ':', field_name_field, '/', field_name_rho
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="rho_bar" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', field_name_initial_conditions_rho_bar
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="rho_bar_z" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', field_name_initial_conditions_rho_bar_z
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="ubc" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', field_name_initial_conditions_ubc
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="dubcdz" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', field_name_initial_conditions_dubcdz
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '            </Grid>'

end subroutine write_xdmf_output_timestep

subroutine write_xdmf_restart_file( xdmf_file_name, hdf5_file_name )
! Writes an XDMF file that represents the contents of the checkpoint restart
! HDF5 file written by the solver.  The file generated, along with the field
! file, may be visualized by anything that understands XDMF data including
! Paraview, VisIt, Mayavi, anything VTK-based, etc.  Creates field variable
! references for each of the internal state variables.
!
! NOTE: The XML written is not validated before it is generated and may
!       be invalid if care is not taken.
!
! NOTE: We jump through hoops to describe the output data as 2D when the
!       number of transverse planes is 1.  This is done to avoid problems with
!       VTK-based visualization tools (specifically, Paraview 5.1.2) which
!       crash when a 3D volume's leading dimension is 1.  Should this not be a
!       VTK issue and rather Greg's inability to understand the intracies of
!       XDMF v2, VTK, and ParaView then the mess below can be greatly
!       simplified.
!
! NOTE: We are not reporting the constant field arrays rho_bar, ubc and
!       dubcdz. Consider modifying this code in the future to add them.
!       Then again, unlike the velocity, we don't really need to visualize
!       the constant fields, rather make sure that they are present in the file.
!
! May 2017
! Gustavo Rivera with Greg Thomsen
!

  use constants, only:  n, nsubx, nsuby, nsubz
  use io_restart, only: restart_name_dubcdz, &
                        restart_name_rho0, restart_name_rho1, restart_name_rho2, &
                        restart_name_rho_bar, restart_name_rho_bar_z, restart_name_ubc, &
                        restart_name_ux0, restart_name_ux1, restart_name_ux2, &
                        restart_name_uy0, restart_name_uy1, restart_name_uy2, &
                        restart_name_uz0, restart_name_uz1, restart_name_uz2

  implicit none

  character(len=*), intent(in) :: xdmf_file_name
  character(len=*), intent(in) :: hdf5_file_name

  ! Dimensions of our data as found in the restart HDF5 file.
  character(len=32)            :: dimensions_string

  integer, parameter           :: unit = 1024

  call get_xdmf_dimensions_string( nsubx * n, nsuby, nsubz * n, dimensions_string )

  open( unit, file=xdmf_file_name )

  call write_xdmf_header( unit )

  ! Write out the grid and its topology.
  call write_xdmf_file_geometry( unit, hdf5_file_name )

  write( unit, '(A)' )              ''
  write( unit, '(A)' )              '        <Grid Name="restart" CollectionType="Spatial" GridType="Collection">'
  write( unit, '(A)' )              '            <Grid Name="solver_conditions">'
  write( unit, '(A)' )              ''
  write( unit, '(A)' )              '                <xi:include xpointer="element(/1/1/1)"/>'
  write( unit, '(A)' )              '                <xi:include xpointer="element(/1/1/2)"/>'
  write( unit, '(A)' )              ''
  write( unit, '(A)' )              '                <Attribute Name="ux0" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_ux0
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="ux1" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_ux1
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="ux2" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_ux2
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  if ( nsuby > 1 ) then

     write( unit, '(A)' )           '                <Attribute Name="uy0" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )       '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )     '                    ', hdf5_file_name, ':', restart_name_uy0
     write( unit, '(A)' )           '                    </DataItem>'
     write( unit, '(A)' )           '                </Attribute>'
     write( unit, '(A)' )           '                <Attribute Name="uy1" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )       '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )     '                    ', hdf5_file_name, ':', restart_name_uy1
     write( unit, '(A)' )           '                    </DataItem>'
     write( unit, '(A)' )           '                </Attribute>'
     write( unit, '(A)' )           '                <Attribute Name="uy2" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )       '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )     '                    ', hdf5_file_name, ':', restart_name_uy2
     write( unit, '(A)' )           '                    </DataItem>'
     write( unit, '(A)' )           '                </Attribute>'
  end if

  write( unit, '(A)' )              '                <Attribute Name="uz0" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_uz0
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="uz1" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_uz1
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="uz2" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_uz2
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="rho0" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_rho0
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="rho1" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_rho1
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="rho2" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_rho2
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="rho_bar" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_rho_bar
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="rho_bar_z" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_rho_bar_z
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  ! XXX: Uncomment the background current information when it is written as
  !      proper 3D data.
  write( unit, '(A)' )              '<!--'
  write( unit, '(A)' )              '                <Attribute Name="ubc" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_ubc
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="dubcdz" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', restart_name_dubcdz
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '-->'

  write( unit, '(A)' )              '            </Grid>'

  call write_xdmf_footer( unit )

  close( unit )

end subroutine write_xdmf_restart_file

subroutine get_xdmf_dimensions_string( nsubx, nsuby, nsubz, dimensions_string )
! Formats the grid's dimensions into a whitespace delimited string.  The
! caller may request that dimensions only be provided for a single transverse
! plane instead of all transverse planes.

  implicit none

  integer, intent(in)              :: nsubx
  integer, intent(in)              :: nsuby
  integer, intent(in)              :: nsubz

  ! Sized for ~1e30 elements. 10-characters per dimension is almost 1e10 each,
  ! and two spaces separating them.
  character(len=32), intent(inout) :: dimensions_string

  ! Create a format specifier string that has exactly the correct number of
  ! digits needed to format our XML output.  This is the form of either
  ! '(INx,A,INz)' for 2D or '(INy,A,INx,A,INz)' for 3D.
  if ( nsuby == 1 ) then
     write( dimensions_string, '(I0,A,I0)' ) nsubx, ' ', nsubz
  else
     write( dimensions_string, '(I0,A,I0,A,I0)' ) nsuby, ' ', nsubx, ' ', nsubz
  endif

end subroutine get_xdmf_dimensions_string

subroutine write_xdmf_header( unit )
! Writes out the "header" for an version 2 XDMF file.  This provides the XML
! declaration (along with a document type and DTD) and opens the XDMF document
! and main computational domain's tags.

  implicit none

  integer, intent(in) :: unit

  write( unit, '(A)' ) '<?xml version="1.0" ?>'
  write( unit, '(A)' ) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'

  ! Create a temporal collection that describes our field outputs.
  write( unit, '(A)' ) '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  write( unit, '(A)' ) '    <Domain>'

end subroutine write_xdmf_header

subroutine write_xdmf_footer( unit )
! Writes out the "footer" for an XDMF file.  This closes the main grid, the
! main computational domain, and the XDMF document itself.

  implicit none

  integer, intent(in) :: unit

  ! Wrap things up.
  write( unit, '(A)' ) '        </Grid>'
  write( unit, '(A)' ) '    </Domain>'
  write( unit, '(A)' ) '</Xdmf>'

end subroutine write_xdmf_footer

subroutine write_xdmf_file_geometry( unit, hdf5_file_name )
! Writes out the XDMF topology and grid for the field file.  The caller can
! request a 2D planar geometry instead of the full 3D geometry.

  use constants, only: n, nsubx, nsuby, nsubz
  use io_field, only:  field_name_x, field_name_y, field_name_z

  implicit none

  integer, intent(in)          :: unit
  character(len=*), intent(in) :: hdf5_file_name

  character(len=32)            :: dimensions_string

  call get_xdmf_dimensions_string( nsubx * n, nsuby, nsubz * n, dimensions_string )

  ! Write out the grid's geometry.
  !
  ! NOTE: Be careful about wrapping lines at 132 columns.
  if ( nsuby == 1 ) then
     write( unit, '(A)' )      '        <!--'
     write( unit, '(A)' )      '          Individual node locations are stored slowest to fastest (X then Z), while the'
     write( unit, '(A)' )      '          order of individual coordinates match the solver''s view of the world (the X-Z'
     write( unit, '(A)' )      '          plane maps to the X-Y plane).'
     write( unit, '(A)' )      '        -->'
     write( unit, '(A,A,A)' )  '        <Topology TopologyType="2DSMesh" NumberOfElements="', &
                                        trim( dimensions_string ), '"/>'
     write( unit, '(A)' )      '        <Geometry GeometryType="X_Y">'
  else
     write( unit, '(A)' )      '        <!--'
     write( unit, '(A)' )      '          Individual node locations are stored slowest to fastest (Y, X, then Z,'
     write( unit, '(A)' )      '          AKA by planes), while the order of individual coordinates match the solver''s'
     write( unit, '(A)' )      '          view of the world (X-Z planes map to X-Y planes).'
     write( unit, '(A)' )      '        -->'
     write( unit, '(A,A,A)' )  '        <Topology TopologyType="3DSMesh" NumberOfElements="', &
                                        trim( dimensions_string ), '"/>'
     write( unit, '(A)' )      '        <Geometry GeometryType="X_Y_Z">'
  end if

  write( unit, '(A,A,A)' )     '            <DataItem Dimensions="', trim( dimensions_string ), &
                                                                     '" NumberType="Float" Precision="8" Format="HDF">'
  write( unit, '(A,A,A,A)' )   '            ', hdf5_file_name, ':', field_name_x
  write( unit, '(A)' )         '            </DataItem>'
  write( unit, '(A,A,A)' )     '            <DataItem Dimensions="', trim( dimensions_string ), &
                                                                     '" NumberType="Float" Precision="8" Format="HDF">'
  write( unit, '(A,A,A,A)' )   '            ', hdf5_file_name, ':', field_name_z
  write( unit, '(A)' )         '            </DataItem>'
  if ( nsuby > 1 ) then
     write( unit, '(A,A,A)' )  '            <DataItem Dimensions="', trim( dimensions_string ), &
                                                                       '" NumberType="Float" Precision="8" Format="HDF">'
     write( unit, '(A,A,A,A)' )'            ', hdf5_file_name, ':', field_name_y
     write( unit, '(A)' )      '            </DataItem>'
  end if
  write( unit, '(A)' )         '        </Geometry>'

end subroutine write_xdmf_file_geometry
