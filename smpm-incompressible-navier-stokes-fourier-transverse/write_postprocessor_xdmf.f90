subroutine write_xdmf_postprocessor_file( xdmf_file_name, hdf5_file_name)
! Writes an XDMF file that represents the contents of the HDF5 postprocessed
! field file written by the solver.  The file generated
! may be visualized by anything that understands XDMF data including
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
! May 2017
! Gustavo Rivera with Greg Thomsen

  implicit none

  character(len=*), intent(in)        :: xdmf_file_name
  character(len=*), intent(in)        :: hdf5_file_name

  integer, parameter                  :: unit = 1024

  open( unit, file=xdmf_file_name )

  call write_xdmf_header( unit )

  ! Write out the grid and its topology.
  call write_xdmf_postprocessorfile_geometry( unit, hdf5_file_name )

  ! Write out all the fields in the postprocessor file.
  call write_xdmf_postprocessor_field( unit, hdf5_file_name )

  call write_xdmf_footer( unit )

  close( unit )

end subroutine write_xdmf_postprocessor_file

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

subroutine write_xdmf_postprocessor_field( unit, hdf5_file_name )
! Writes out field variables for a single timestep in XDMF.  This correctly
! handles writing out the initial conditions as step #0.

  use constants, only: n, nsubz, nsuby, nsubx
  use io_post, only:   post_name_div, &
                       post_name_drhodx, post_name_drhody, post_name_drhodz, &
                       post_name_dudx, post_name_dudy, post_name_dudz, &
                       post_name_dvdx, post_name_dvdy, post_name_dvdz, &
                       post_name_dwdx, post_name_dwdy, post_name_dwdz, &
                       post_name_enstrophy, post_name_kinetic_energy, &
                       post_name_omega_x, post_name_omega_y, post_name_omega_z, &
                       post_name_p, post_name_rho, post_name_rho_bar, &
                       post_name_stream, &
                       post_name_ux, post_name_uy, post_name_uz

  implicit none

  integer, intent(in)          :: unit
  character(len=*), intent(in) :: hdf5_file_name

  character(len=32)            :: dimensions_string

  call get_xdmf_dimensions_string( nsubx * n, nsuby, nsubz * n, dimensions_string )

  ! Write out all the postprocessed field.

  write( unit, '(A)' )              ''
  write( unit, '(A)' )              '        <Grid Name="postprocessor" CollectionType="Spatial" GridType="Collection">'
  write( unit, '(A)' )              '            <Grid Name="solver_conditions">'
  write( unit, '(A)' )              ''
  write( unit, '(A)' )              '                <xi:include xpointer="element(/1/1/1)"/>'
  write( unit, '(A)' )              '                <xi:include xpointer="element(/1/1/2)"/>'
  write( unit, '(A)' )              ''
  write( unit, '(A)' )              '                <Attribute Name="div" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_div
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  write( unit, '(A)' )              '                <Attribute Name="drhodx" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_drhodx
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  if ( nsuby > 1) then 
     write( unit, '(A)' )              '                <Attribute Name="drhody" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                         trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_drhody
     write( unit, '(A)' )              '                    </DataItem>'
     write( unit, '(A)' )              '                </Attribute>'
  endif

  write( unit, '(A)' )              '                <Attribute Name="drhodz" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_drhodz
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="dudx" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_dudx
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  if ( nsuby > 1 ) then
     write( unit, '(A)' )              '                <Attribute Name="dudy" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                         trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_dudy
     write( unit, '(A)' )              '                    </DataItem>'
     write( unit, '(A)' )              '                </Attribute>'
  endif

  write( unit, '(A)' )              '                <Attribute Name="dudz" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_dudz
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  
  if ( nsuby > 1 ) then
     write( unit, '(A)' )              '                <Attribute Name="dvdx" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                         trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_dvdx
     write( unit, '(A)' )              '                    </DataItem>'
     write( unit, '(A)' )              '                </Attribute>'
     write( unit, '(A)' )              '                <Attribute Name="dvdy" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                         trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_dvdy
     write( unit, '(A)' )              '                    </DataItem>'
     write( unit, '(A)' )              '                </Attribute>'
     write( unit, '(A)' )              '                <Attribute Name="dvdz" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                         trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_dvdz
     write( unit, '(A)' )              '                    </DataItem>'
     write( unit, '(A)' )              '                </Attribute>'
  endif

  write( unit, '(A)' )              '                <Attribute Name="dwdx" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_dwdx
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  if ( nsuby > 1 ) then
     write( unit, '(A)' )              '                <Attribute Name="dwdy" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                         trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_dwdy
     write( unit, '(A)' )              '                    </DataItem>'
     write( unit, '(A)' )              '                </Attribute>'
  endif

  write( unit, '(A)' )              '                <Attribute Name="dwdz" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_dwdz
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'
  
  write( unit, '(A)' )              '                <Attribute Name="enstrophy" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_enstrophy
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="KE" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_kinetic_energy
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  if ( nsuby > 1 ) then
     write( unit, '(A)' )              '                <Attribute Name="omega_x" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                         trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_omega_x
     write( unit, '(A)' )              '                    </DataItem>'
     write( unit, '(A)' )              '                </Attribute>'
  endif

  write( unit, '(A)' )              '                <Attribute Name="omega_y" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_omega_y
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  if ( nsuby > 1 ) then
     write( unit, '(A)' )              '                <Attribute Name="omega_z" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                         trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_omega_z
     write( unit, '(A)' )              '                    </DataItem>'
     write( unit, '(A)' )              '                </Attribute>'
  endif

  write( unit, '(A)' )              '                <Attribute Name="pressure" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_p
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="rho" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_rho
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="rho_bar" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_rho_bar
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="ux" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_ux
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  if ( nsuby > 1 ) then
     write( unit, '(A)' )              '                <Attribute Name="uy" AttributeType="Scalar" Center="Node">'
     write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                         trim( dimensions_string ), '" Format="HDF">'
     write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_uy
     write( unit, '(A)' )              '                    </DataItem>'
     write( unit, '(A)' )              '                </Attribute>'
  endif

  write( unit, '(A)' )              '                <Attribute Name="uz" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_uz
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '                <Attribute Name="streamfunction" AttributeType="Scalar" Center="Node">'
  write( unit, '(A,A,A)' )          '                    <DataItem DataType="Float" Precision="8" Dimensions="', &
                                                                      trim( dimensions_string ), '" Format="HDF">'
  write( unit, '(A,A,A,A)' )        '                    ', hdf5_file_name, ':', post_name_stream
  write( unit, '(A)' )              '                    </DataItem>'
  write( unit, '(A)' )              '                </Attribute>'

  write( unit, '(A)' )              '            </Grid>'

end subroutine write_xdmf_postprocessor_field



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

subroutine write_xdmf_postprocessorfile_geometry( unit, hdf5_file_name )
! Writes out the XDMF topology and grid for the post file.  The caller can
! request a 2D planar geometry instead of the full 3D geometry.

  use constants,  only:        n, nsubx, nsuby, nsubz
  use io_post, only:           post_name_x, post_name_y, post_name_z

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
  write( unit, '(A,A,A,A)' )   '            ', hdf5_file_name, ':', post_name_x
  write( unit, '(A)' )         '            </DataItem>'
  write( unit, '(A,A,A)' )     '            <DataItem Dimensions="', trim( dimensions_string ), &
                                                                     '" NumberType="Float" Precision="8" Format="HDF">'
  write( unit, '(A,A,A,A)' )   '            ', hdf5_file_name, ':', post_name_z
  write( unit, '(A)' )         '            </DataItem>'
  if ( nsuby > 1 ) then
     write( unit, '(A,A,A)' )  '            <DataItem Dimensions="', trim( dimensions_string ), &
                                                                       '" NumberType="Float" Precision="8" Format="HDF">'
     write( unit, '(A,A,A,A)' )'            ', hdf5_file_name, ':', post_name_y
     write( unit, '(A)' )      '            </DataItem>'
  end if
  write( unit, '(A)' )         '        </Geometry>'

end subroutine write_xdmf_postprocessorfile_geometry
