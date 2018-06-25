#!/bin/sh

# script that generates a list of intermediate files that can be ignored by
# Mercurial from Fortran source code found in a user-specified directory.

# display the usage on standard output.
print_usage()
{
    echo "Usage: $0 [-E <extension>] [-h] [<directory>]"
    echo
    echo "Mercurial ignore file generation script.  Generates a list of object and module"
    echo "file names from the Fortran source code located in <directory>, or the current"
    echo "directory if it is omitted.  The list of files to ignore is displayed on standard"
    echo "output and is suitable for inclusion in a Mercurial .hgignore file."
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

# by default, work in the current directory.
DIRECTORY=.

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

# figure out how the user called us.
if [ $# -gt 1 ]; then
    echo "Too many arguments provided.  Expected 0 or 1, received $#."
    print_usage
    exit 1
elif [ $# -eq 1 ]; then
    DIRECTORY=$1
fi

# ensure that the directory exists, otherwise we have nothing to do.
if [ ! -d "${DIRECTORY}" ]; then
    echo "${DIRECTORY} doesn't exist.  Nothing to do."
    exit 1
fi

# Fortran source code compiles into object files that should be ignored.
echo "# ignore object files that correspond to known source files."

for F90 in ${DIRECTORY}/*.${EXTENSION}; do
    # skip over ifort's generated interface files.
    if /bin/echo ${F90} | grep -q "__genmod.${EXTENSION}$"; then
        continue
    fi

    echo ${F90} | sed -e "s/.${EXTENSION}\$/.o/"
done | sort | sed -e "s#^${DIRECTORY}/##"

# Fortran modules also create .mod files that should be ignored as well.
echo
echo "# ignore module files that correspond to known module interfaces."

for F90 in ${DIRECTORY}/mod_*.${EXTENSION}; do
    echo ${F90} | sed -e 's/mod_//' -e "s/.${EXTENSION}/.mod/"
done | sort | sed -e "s#^${DIRECTORY}/##"

# handle Intel's ifort generates interface files separately from normal
# Fortran source, since interface files are not generated for module files.
echo
echo "# ignore Intel ifort's interface files."

# we make two passes through the files to create a list of interfaces to
# ignore.  the first is based off of the source file names, which covers
# functions and sub-routines.  the second is based off of sub-routines
# found within each source file.  we generate everything in a subshell
# so that we don't have to worry about duplicating ignore rules within
# each part.
(for F90 in ${DIRECTORY}/*.${EXTENSION}; do
    # filter out a few things that shouldn't have ignore rules generated
    # for them:
    #
    #   1. module files
    #   2. interface files
    #   3. main.${EXTENSION} since it's likely to contain a subprogram
    #      definition
    if /bin/echo ${F90} | egrep -q -e "(/mod_|(main|__genmod).${EXTENSION})"; then
        continue
    fi

    echo ${F90} | sed -e "s/.${EXTENSION}/__genmod.${EXTENSION}/"
    echo ${F90} | sed -e "s/.${EXTENSION}/__genmod.mod/"
done

# sub-routines in each of the source files (but not the auto-generated
# ones) also get interface files.
#
# NOTE: we suppress egrep's output prefix (-h) so that we can properly handle
#       cases where the subroutine closure statement does not start in the
#       first column.  without this, the first egrep could produce things
#       like:
#
#         file.f90:     end subroutine some_subroutine
#
#       which causes the awk to pick out the wrong sub-routine name.
for SUBROUTINE in `egrep -hi -e 'end.*subroutine' ${DIRECTORY}/*.${EXTENSION} | \
                     grep -v __genmod | awk '{ print tolower($3) }'`; do
    echo "${SUBROUTINE}__genmod.${EXTENSION}"
    echo "${SUBROUTINE}__genmod.mod"
done

# remove the directory before sorting since the subroutine search does not
# have the directory prefix.  otherwise, we get spurious duplicates.
) | sed -e "s#^${DIRECTORY}/##" | sort -u
