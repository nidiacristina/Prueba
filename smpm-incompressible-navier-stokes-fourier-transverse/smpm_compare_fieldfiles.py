#!/usr/bin/env python3

# takes two HDF5 field files generated from the SMPM solver and compares
# thems to determine if they solved the same problem to within tolerance.
# this is intended to validate different solver implementations/configurations
# against known good outputs.

import getopt
import subprocess
import sys

# option strings used to specify which portions of the field file to use in
# the comparison.
OPTION_COMPARE_CONFIG = "config"
OPTION_COMPARE_OUTPUT = "output"

# our solver operates in double precision so we set a fairly tight relative
# tolerance by default.
DEFAULT_TOLERANCE = 1e-10

HDF5_DIFF = "h5diff"

def usage( script_name ):
    """
    Takes a name of the script (full path, name, etc) and prints its usage to standard
    output.
    """

    print( """Usage: {script_name:s} [-d] [-h] [-q <content>[,<content>[...]] [-v] <file1> <file2> [<tolerance>]

Takes two HDF5 files produced by the SMPM Incompressible Navier-Stokes solver
and compares their contents to determine if they have solved the same problem
to with a specified tolerance.  This allows validation of different solver
implementations/configurations against a known good output.

The contents of <file2> are compared against those of <file1> according to the
<content> requested.  By default, both the problem configuration and the field
output are compared though either may be compared individually.  The HDF5 groups
and datasets compared are as follows for <content>:

  {configuration_option:s}

    /grid/            - Grid parameters
    /configuration/   - Solver configuration (density, viscosity, facrobin, ...)
    /field/step0/     - Initial conditions

   NOTE: Only the datasets related to the problem's physics in the groups above
         are compared.  Execution-dependent datasets (number of processors, wall-clock
         time, etc) are ignored.

  {field_option:s}

    /field/stepN/     - Field outputs.

The command line options above are described below:

  -d                            Debug the comparison process.  Prints the h5diff
                                commands before they're executed.
  -h                            Displays this help message and exits.
  -q <content>[,<content>[...]  Specifies which content to compare.  See above for
                                the datasets compared.  Must be either
                                '{configuration_option:s}' or '{field_option:s}'.
                                If omitted, defaults to '{configuration_option:s},{field_option:s}'.
  -v                            Be verbose and list the datasets that differ and
                                how many elements within each.
""".format( script_name=script_name,
            configuration_option=OPTION_COMPARE_CONFIG,
            field_option=OPTION_COMPARE_OUTPUT ) )

def build_compare_command( file_name_1, file_name_2, tolerance,
                           dataset, excluded_paths=[] ):
    """
    Builds a command line list to compare a common dataset between HDF5 files.
    Specific datasets may be omitted when comparing a common group.

    Takes 5 arguments:

      file_name_1   - First HDF5 file name.
      file_name_2   - Second HDF5 file name.
      tolerance     - Relative tolerance to compare the common dataset with.
      dataset       - Path to common dataset within both HDF5 files.  May be
                      absolute or relative.
      excluded_path - Path to datasets within dataset to exclude during
                      comparison.  This is a relative path.

    Returns 1 value:

      command_list -
    """

    import itertools

    # NOTE: be careful when specifying paths to h5diff.  if a group has a trailing
    #       slash (/) appended, then any paths generated from it *MUST* not have
    #       a leading slash (/) prepended otherwise the target dataset/group will
    #       not be handled properly.
    #
    #       for example, if the collocation dimension (/grid/n) should be ignored
    #       then the following paths will properly refer to it:
    #
    #           /grid    /n    ->   /grid/n
    #           /grid/    n    ->   /grid/n
    #
    #       and the following will *not* refer to it:
    #
    #           /grid/   /n    ->   /grid//n
    #
    #       currently we don't attempt to fix this for the user, we simply
    #       warn them and let them deal with the issue.
    if len( excluded_paths ) > 0:
        if ((dataset.endswith( "/" ) and sum( map( lambda x: x.startswith( "/" ), excluded_paths ) ) > 0) or
            (not dataset.endswith( "/" ) and
             sum( map( lambda x: x.startswith( "/" ), excluded_paths ) ) != len( excluded_paths ))):
            print( ("Requested dataset ({dataset:s}) and the exclusion " +
                    "paths ({exclusion:s}) are incompatible.  Please review!").format( dataset=dataset,
                                                                                       exclusion=", ".join( excluded_paths ) ),
                   file=sys.stderr )

    # base command that we build from.  this requests comparing the entirety of
    # both files' contents until we ask for specific datasets/groups instead,
    # which we do for each of the executions of h5diff below.
    command = [HDF5_DIFF]
    command.extend( itertools.chain.from_iterable( map( lambda x: ("--exclude-path", x), excluded_paths ) ) )
    command.extend( ["-p", "{:.15f}".format( tolerance ),
                     file_name_1,
                     file_name_2,
                     dataset] )

    return command

def compare_grid( file_name_1, file_name_2, tolerance, debug_flag ):
    """
    Compares the grid of each file to one another to determine if they're the
    same.

    Takes 3 arguments:

      file_name_1 - First SMPM field file name.
      file_name_2 - Second SMPM field file name.
      tolerance   - Tolerance to verify similarity with.

    Returns 1 value:

      similar_flag - Boolean indicating whether the files' grids are similar.
    """

    dimensions_command = build_compare_command( file_name_1,
                                                file_name_2,
                                                tolerance,
                                                "/grid",
                                                excluded_paths=("/x", "/y", "/z") )

    grid_command = build_compare_command( file_name_1,
                                          file_name_2,
                                          tolerance,
                                          "/grid",
                                          excluded_paths=("/n", "/mx", "/my", "/mz") )

    if not compare_and_summarize( dimensions_command, debug_flag ):
        return False
    return compare_and_summarize( grid_command, debug_flag )

def compare_configuration( file_name_1, file_name_2, tolerance, debug_flag ):
    """
    Compares the configuration of each file to one another to determine if
    they're the same.

    Takes 3 arguments:

      file_name_1 - First SMPM field file name.
      file_name_2 - Second SMPM field file name.
      tolerance   - Tolerance to verify similarity with.

    Returns 1 value:

      similar_flag - Boolean indicating whether the files' configurations are
                     similar.
    """

    # compare the configurations, sans GMRES parameters.  we skip GMRES
    # thresholds since they could be set high enough, yet different from
    # each other, where neither invocation exceeded them, yet converged to
    # the appropriate tolerance.
    #
    # we also ignore the kernel vectors for the capacitance matrix of the
    # Poisson operator (uC) and the Poisson operator (uL) itself.  these can
    # span the same space but be completely different run to run so we don't
    # bother comparing them.
    # XXX: is this true?
    config_command = build_compare_command( file_name_1,
                                            file_name_2,
                                            tolerance,
                                            "/configuration",
                                            excluded_paths=("/gmres", "/uC", "/uL") )

    # compare the initial conditions in their entirety.
    ic_command = build_compare_command( file_name_1,
                                        file_name_2,
                                        tolerance,
                                        "/field/step0" )

    if not compare_and_summarize( config_command, debug_flag ):
        return False
    return compare_and_summarize( ic_command, debug_flag )

def compare_field( file_name_1, file_name_2, tolerance, debug_flag ):
    """
    Compares the field outputs of each file to one another to determine if
    they're the same.

    Takes 3 arguments:

      file_name_1 - First SMPM field file name.
      file_name_2 - Second SMPM field file name.
      tolerance   - Tolerance to verify similarity with.

    Returns 1 value:

      similar_flag - Boolean indicating whether the files' configurations are
                     similar.
    """

    # compare the entire field while taking care to ignore the initial
    # conditions (step0).
    field_command = build_compare_command( file_name_1,
                                           file_name_2,
                                           tolerance,
                                           "/field",
                                           excluded_paths=("/step0",) )
    return compare_and_summarize( field_command, debug_flag )

def compare_and_summarize( command, debug_flag ):
    """
    Executes a comparison command and summarizes the output if there are
    differences.  Returns a flag indicating whether the comparison resulted in
    similar files or dissimilar files.  If the comparison does not result in
    similar files, either a summary of differences or the error that occurred
    is printed to standard output.

    Takes 1 arguments:

      command - Iterable containing the comparison command's components.
                Must be suitable to pass to subprocess.check_output().

    Returns 1 value:

      similar_flag - Boolean indicating whether the comparison command
                     resulted in similarity (True) or dissimilarity (False).
    """

    if debug_flag:
        print( " ".join( command ) )

    try:
        output = subprocess.check_output( command )
    except subprocess.CalledProcessError as e:
        if e.returncode == 1:
            print( parse_h5diff_output( e.output ) )
        else:
            print( "Error code {:d}.".format( e.returncode ) )

        return False

    return True

def parse_h5diff_output( output_bytes ):
    """
    Parses the output of a non-verbose h5diff invocation and summarizes the output
    under the assumption that we're comparing a common dataset between two
    different files.  In this case the two dataset names are identical and we
    only need how many differences were found.  Information not related to
    differences (e.g. indicators of non-comparable datasets) are ignored.

    Takes 1 argument:

      output_bytes - UTF-8 encoded byte stream from an h5diff comparison.

    Returns 1 value:

      summary - String containing the simplified summary.

    """

    # name of the most recently compared dataset.
    dataset_name = ""

    # walk through the output looking for output indicating comparable datasets
    # that had differences.  these have the form of:
    #
    #   dataset: <foo> and <bar>
    #   N differences found.
    #
    # we break this down to track <foo> and N so a simplified summary can be
    # generated.  everything else will be ignored.
    summary = ""
    for line in output_bytes.decode( "utf-8" ).split( "\n" )[:-1]:
        words = line.split( " " )
        if line.startswith( "dataset:" ):
            # found the start of a new difference block, of the form:
            #
            #    dataset: <foo> and <bar>
            #
            # track the name of the dataset in the first file compared.
            dataset_name = line.split( " " )[1][1:-1]
        elif words[1:] == ["differences", "found"]:
            # found the end of a difference block, of the form:
            #
            #    N differences found.
            #
            # use the count to update our summary output.
            count = int( words[0] )
            summary += "{:s} had {:d} difference{:s}.\n".format( dataset_name,
                                                                 count,
                                                                 "" if count == 0 else "s" )

    return summary

def main( argv ):
    """
    Takes a list of arguments, including the name of the application calling this method,
    and compares two field files for correctness within a given tolerance.
    """

    # by default we compare both configuration and field outputs.  since users
    # can only request a comparison, we set both to false and detect a lack of
    # user request later.
    compare_configuration_flag = False
    compare_output_flag        = False

    # we assume normal operations unless requested by the user.
    debug_flag = False

    # parse our command line options.
    try:
        opts, args = getopt.getopt( argv[1:], "dhq:v" )
    except getopt.GetoptError as error:
        print( "Error processing option: {:s}\n".format( str( error ) ),
               file=sys.stderr )
        sys.exit( 1 )

    # handle any valid options were were presented.
    for opt, arg in opts:
        if opt == '-d':
            debug_flag = True
        if opt == '-h':
            usage( argv[0] )
            sys.exit()
        if opt == '-q':
            # validate that the content type specified.
            options = (OPTION_COMPARE_CONFIG, OPTION_COMPARE_OUTPUT)

            if not (arg in options):
                print( ("Unknown content type specified ({specified:s}).  " +
                        "Must be one of {content:s}.\n").format( specified="'" + arg + "'",
                                                                 content=" or ".join( map( lambda x: "'" + x + "'",
                                                                                           options ) ) ) )
                sys.exit( 1 )
            elif arg == OPTION_COMPARE_CONFIG:
                compare_configuration_flag = True
            elif arg == OPTION_COMPARE_OUTPUT:
                compare_output_flag = True

    # make sure we were called properly.
    if not (1 < len( args ) < 4):
        print( "Expected 2 or 3 arguments but received {:d}.\n".format( len( args ) ),
               file=sys.stderr )
        usage( argv[0] )
        sys.exit( 1 )

    file_name_1 = args[0]
    file_name_2 = args[1]
    if len( args ) == 3:
        tolerance = float( args[2] )
    else:
        tolerance = DEFAULT_TOLERANCE

    # we compare everything if the user has not requested a specific type of
    # content compared.
    if not (compare_configuration_flag or compare_output_flag):
        compare_configuration_flag = True
        compare_output_flag        = True

    # regardless of user request, we always compare the grids first.  we can
    # save ourselves a lot of work if we identify different grids up front as
    # we can't possibly have the same field outputs.
    if not compare_grid( file_name_1, file_name_2, tolerance, debug_flag ):
        print( "{file1:s} and {file2:s} have different grids.\n".format( file1=file_name_1,
                                                                         file2=file_name_2 ),
               file=sys.stderr )
        sys.exit( 1 )

    if compare_configuration_flag:
        result = compare_configuration( file_name_1, file_name_2, tolerance, debug_flag )

        if not result:
            print( "{file1:s} and {file2:s} appear to have different configurations.\n".format( file1=file_name_1,
                                                                                                file2=file_name_2 ),
                   file=sys.stderr )
            sys.exit( 1 )

    if compare_output_flag:
        result = compare_field( file_name_1, file_name_2, tolerance, debug_flag )

        if not result:
            print( "{file1:s} and {file2:s} appear to have different outputs.\n".format( file1=file_name_1,
                                                                                         file2=file_name_2 ),
                   file=sys.stderr )
            sys.exit( 1 )

    # we compared everything requested and found no differences.
    sys.exit( 0 )

if __name__ == "__main__":
    main( sys.argv )
