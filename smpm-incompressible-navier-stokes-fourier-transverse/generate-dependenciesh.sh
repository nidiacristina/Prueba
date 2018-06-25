#!/bin/sh

# first order attempt at generating dependency information for Fortran code
# suitable for use within a Makefile.

# display the usage on standard output.
print_usage()
{
    echo "Usage: $0 [-E <extension>] [-h] <filename> [<filename> [...]]"
    echo
    echo "Fortran dependency generation script.  Takes a list of Fortran source code file"
    echo "names and generates dependency information suitable for inclusion within a"
    echo "Makefile on standard output."
    echo
    echo "The command line flags shown above are described below:"
    echo
    echo "    -E <extension>  Process files that have <extension>.  If omitted, defaults"
    echo "                    to \"${EXTENSION}\"."
    echo "    -h              Display this help message and exit."
    echo
}

# by default, operate on Fortran 90 files.
EXTENSION=f90

while getopts "E:h" FLAGS
do
    case ${FLAGS} in
        E)
            EXTENSION=${OPTARG}
            ;;
        h)
            print_usage
            exit 0
            ;;
    esac
done

shift `expr ${OPTIND} - 1`

if [ $# -eq 0 ]; then
    echo "No source files provided!  Exiting."
    exit 1
fi

# walk through each source file and generate dependency information.  note that
# regardless of the order provided on the command line, all of the output is
# sorted.  care is also taken to ensure that no trailing whitespace is present
# which occurs when a file does not have any dependencies.
for SOURCE in $*; do

    # start by mapping the source file to the equivalent object file, though
    # do not terminate the line with a new line.
    /bin/echo -n ${SOURCE} | sed -e "s/.${EXTENSION}$/.o: /";

    # identify each of the Fortran module inclusion lines, ignoring certain
    # official modules that aren't ours (OpenMP, HDF5, MPI, Fortran 2003 C
    # bindings, etc), map each of them to the corresponding module's object
    # file name.  all dependencies are sorted uniquely (to avoid duplication
    # when files contain multiple subroutines) and combined onto a single line
    # so as to make it easy to see changes in the future.
    #
    # NOTE: the initial egrep's regular expression is not precise enough
    #       and occasionally matching things like "because f".  to compensate,
    #       we assume that these would be in comment strings, and simply
    #       filter out lines that have comments.  so far, this has worked
    #       without problem.
    egrep -e '[[:space:]]*use [^[[:space:]]]*' ${SOURCE} | grep -v '\!' | \
        egrep -v -e '(omp_lib|HDF5|mpi|iso_c_binding)' | \
        sed -e 's/^.*use //' -e 's/ *$//' -e 's/, *only.*$//' -e 's/^/mod_/' -e 's/$/.o/' | \
        sort -u | \
        xargs ;
done | sed -e 's/ *$//' |
       sort
