# collection of shell functions and  useful for executing the SMPM Incompressible
# Navier-Stokes solver.

# actions available for testing frameworks, intended to reduce the dependencies
# necessary on the execution environment which typically does not provide all of
# the tools needed for preparing SMPM simulation environments.
#
# preparation does everything needed for execution though does not invoke a
# solver or test driver.  execution assumes that preparation has already been
# completed and only invokes the solver or test driver.
WOODBURY_ACTION_PREPARE="prepare"
WOODBURY_ACTION_EXECUTE="execute"

# takes a directory, searches it for valid test cases immediately below it,
# and prints out the relative directories of the unique cases found.
enumerate_available_test_cases()
{
    SEARCH_DIRECTORY="$1"

    # if the directory provided does not end a trailing slash ('/'), then append
    # one.  it is needed so that the sed substitution below properly strips leading
    # directories.
    if echo "${SEARCH_DIRECTORY}" | grep -v -q '/$' ; then
        SEARCH_DIRECTORY="${SEARCH_DIRECTORY}/"
    fi

    # enumerate the directories containing input file generators.  we only
    # look for them in a directory within the search directory.  symlinks are
    # followed so we have flexibility in constructing our test cases.
    find -L "${SEARCH_DIRECTORY}" -mindepth 2 -maxdepth 2 -name '*.m' -type f -exec dirname {} \; | \
        sort -u | \
        sed -e "s#^${SEARCH_DIRECTORY}##"

    unset SEARCH_DIRECTORY
}

# takes a directory, searches it for a file (default-tests.txt) containing a
# list of default test case directories (one per line), and prints out the
# relative directories of the test cases that exist on disk.
enumerate_default_test_cases()
{
    SEARCH_DIRECTORY="$1"

    DEFAULTS_FILE="${SEARCH_DIRECTORY}/default-tests.txt"

    # there are no default test cases if nothing exists in our search
    # directory.
    if [ ! -f "${DEFAULTS_FILE}" ]; then
        return
    fi

    # sanity check the default test cases and only keep the ones that
    # exist on disk.  note that this doesn't verify that this is a
    # valid test case as that is someone else's problem.
    for TEST_CASE in `cat "${DEFAULTS_FILE}"`; do
        if [ -d "${SEARCH_DIRECTORY}/${TEST_CASE}" ]; then
            echo "${TEST_CASE}"
        fi
    done
}

# takes a space delimited list of available test cases, a space delimited list
# of target test cases, and prints out the targets that are available.
# targets that do not correspond to an available test generate a warning.  the
# order of the targets printed is retained relative to the available test case
# list.
select_test_cases()
{
    AVAILABLE_TEST_CASES="$1"
    TARGET_TEST_CASES="$2"

    SELECTED_TEST_CASES=

    for TARGET in ${TARGET_TEST_CASES}; do
        if echo "${AVAILABLE_TEST_CASES}" | egrep -q -e "^${TARGET}$" ; then
            SELECTED_TEST_CASES="${SELECTED_TEST_CASES} ${TARGET}"
        else
            echo "Test case '${TARGET}' is not available.  Skipping!" >&2
        fi
    done

    echo "${SELECTED_TEST_CASES}"

    unset AVAILABLE_TEST_CASES TARGET_TEST_CASES SELECTED_TEST_CASES TARGET
}

# prepares a specified test case for execution by running its input
# generator.  takes a test case's directory, finds the test case generator,
# and executes it.  the global variable INPUT_FILE is set to the generated
# input file on success, or cleared if an error occurred.
prepare_test_case()
{
    DIRECTORY="$1"
    PREPARE_FLAG="$2"

    if [ -n "${PREPARE_FLAG}" -a "${PREPARE_FLAG}" = "yes" ]; then
        # identify the input files generator in this test directory.
        INPUT_GENERATOR=`find -L ${DIRECTORY} -type f -name "*.m"`
        GENERATOR_PREFIX=`echo ${INPUT_GENERATOR} | sed -e 's/\.m$//'`

        # fail this test case if we don't have a single generator.
        GENERATOR_COUNT=`echo "${INPUT_GENERATOR}" | wc -l`
        if [ ${GENERATOR_COUNT} -ne 1 ]; then
            test_failed "Unknown" "Incorrect number of input generators (${GENERATOR_COUNT}) in ${DIRECTORY}."
            return
        fi

        # prepare the case via its generator.
        #
        # NOTE: the number of processors is wonky with respect to both MPI and
        #       OpenMP.
        #
        #       MPI:    an incorrect memory estimate occurs when the number of
        #               subdomains in the test case isn't evenly divisible by the
        #               processor count.
        #       OpenMP: the current OpenMP thread count isn't considered in the
        #               estimate.
        #
        INPUT_FILE=`${MATLAB_COMMAND} "addpath(genpath('${MATLAB_PATH}'));smpm_prepare_test_case('${GENERATOR_PREFIX}',${NUMBER_PROCESSORS});exit;"`
    else
        INPUT_FILE=`find -L ${DIRECTORY} -type f -name "*_in" | head -n1`
    fi

    # if we didn't get output, we didn't prepare things properly and we have a
    # problem.
    if [ -z "${INPUT_FILE}" ]; then
        test_failed "Unknown" "Failed to generate an input file from ${INPUT_GENERATOR}."
    fi
}

# prints the number of cores found on the system.
get_system_core_count()
{
    # currently, detection of the cores only works on Linux.
    if [ `uname` = "Linux" ]; then
        grep MHz /proc/cpuinfo | wc -l
    else
        echo 1
    fi
}

# takes a processor count and outputs "yes" or "no" to standard output if
# the count was valid or not, respectively.
is_valid_processor_count()
{
    COUNT=$1

    if echo "${COUNT}" | egrep -q -e '^[1-9][0-9]*$' ; then
        echo "yes"
    else
        echo "no"
    fi

    unset COUNT
}

# takes a test input file name and echos the number of horizontal subdomains
# to standard output.
get_number_x_subdomains()
{
    INPUT_FILENAME=$1

    grep '^[ ]*nsubx' ${INPUT_FILENAME} | sed -e 's/.*=//'
}

# takes a core count and a test input file name and echos the maximum number
# of cores that can be used to successfully execute the solver on the input.
# note that the number of cores returned may be fewer than the number
# specified.
get_suitable_core_count()
{
    NUMBER_CORES=$1
    INPUT_FILENAME=$2

    # pull out the number of horizontal sub-domains from the configuration
    # file.  since these are partitioned and distributed across cores, the
    # number of cores requested from MPI needs to match.
    NUMBER_SUBDOMAINS=`get_number_x_subdomains ${INPUT_FILENAME}`

    # if we cannot figure out how many sub-domains there are, assume a single
    # one.
    if [ -z "${NUMBER_SUBDOMAINS}" ]; then
        echo "Failed to determine the number of sub-domains in ${INPUT_FILENAME}." >&2
        NUMBER_SUBDOMAINS=1
    fi

    # compute the maximum number of cores suitable for the sub-domains in this
    # test.
    compute_suitable_core_count ${NUMBER_CORES} ${NUMBER_SUBDOMAINS}

    unset NUMBER_CORES INPUT_FILENAME NUMBER_SUBDOMAINS
}

# takes a core count and a sub-domain count and computes a suitable core count
# such that each core can be associated with an integral number of
# sub-domains.  the core count computed is echoed to standard out and will be
# no larger than the core count provided.
compute_suitable_core_count()
{
    NUMBER_CORES=$1
    NUMBER_SUBDOMAINS=$2

    NUMBER_SUITABLE_CORES=${NUMBER_CORES}

    # validate that our parameters are sensible before we attempt to find
    # a suitable core count.  no matter what, a core count of 1 is
    # suitable for the solver.
    #
    # NOTE: as of 2014/06/07 the solver does not properly handle a single
    #       core properly.  it is intended that this will be fixed though.
    if [ -z "${NUMBER_SUITABLE_CORES}" -o -z "${NUMBER_SUBDOMAINS}" -o \
         "${NUMBER_SUITABLE_CORES}" -lt 2 ]; then
        NUMBER_SUITABLE_CORES=1
    fi

    # start with the number of cores provided and see if it produces an
    # integral number of sub-domains per core.  if the preferred core count
    # does not produce an integral number of sub-domains per core, it is
    # decremented and checked again.
    while [ ${NUMBER_SUITABLE_CORES} -gt 1 ]; do
        DOMAINS_PER_CORE=`expr ${NUMBER_SUBDOMAINS} / ${NUMBER_SUITABLE_CORES}`

        if [ `expr ${DOMAINS_PER_CORE} '*' ${NUMBER_SUITABLE_CORES}` -eq ${NUMBER_SUBDOMAINS} ] ; then
            break;
        fi

        NUMBER_SUITABLE_CORES=`expr ${NUMBER_SUITABLE_CORES} - 1`
    done

    echo ${NUMBER_SUITABLE_CORES}

    unset NUMBER_CORES NUMBER_SUBDOMAINS NUMBER_SUITABLE_CORES DOMAINS_PER_CORE
}

# takes a variable name and increments the named variable's contents by 1.
increment_variable()
{
    VARIABLE_NAME=$1

    # constructs something like:
    #
    #   FOO=`expr ${FOO} + 1`
    EVAL_STRING="${VARIABLE_NAME}=\`expr \${${VARIABLE_NAME}} + 1\`"

    # have the interpreter evaluate the entire string and update the
    # appropriate variable.
    eval ${EVAL_STRING}
}

# finds a command that can evaluate MATLAB code while respecting the caller's
# wishes to favor Octave or not.  the command suitable for invoking code is
# echoed on standard output.
#
# NOTE: this finds a command suitable for appending a MATLAB command string
#       and successfully executing it, *not* a path to said MATLAB command.
get_matlab_command()
{
    SKIP_MATLAB_FLAG=$1

    # we need to check MATLAB if we're not told to explicitly ignore it.
    if [ "${SKIP_MATLAB_FLAG}" = "no" ]; then
        MATLAB_COMMAND=`which matlab 2>/dev/null`

        # if we found MATLAB, then we're good.
        if [ -n "${MATLAB_COMMAND}" ]; then
            echo "${MATLAB_COMMAND} -r "
            return
        fi
    fi

    # either MATLAB isn't installed or we shouldn't be using it.  see if
    # Octave is available.

    OCTAVE_COMMAND=`which octave 2>/dev/null`

    # if we found Octave, then we're good as well.  make sure we keep Octave
    # from announcing itself and polluting standard output.
    if [ -n "${OCTAVE_COMMAND}" ]; then
        echo "${OCTAVE_COMMAND} --quiet --eval"
    fi

    # whomp, whomp.
    echo ""
}

# parses a list of actions and verifies they are known to the testing framework.
# if they are, this function echoes the actions on standard output, otherwise
# returns an empty string and issues an error on standard error.  if an empty
# action is specified, all of the valid actions are returned.
#
# NOTE: this currently only supports a single action rather than a list.
parse_test_actions()
{
    REQUEST_STRING="${1}"

    # XXX: handle a list of actions properly.
    if [ -z "${REQUEST_STRING}" ]; then
        echo "${WOODBURY_ACTION_PREPARE} ${WOODBURY_ACTION_EXECUTE}"
    elif [ "${REQUEST_STRING}" = "${WOODBURY_ACTION_PREPARE}" ]; then
        echo "${WOODBURY_ACTION_PREPARE}"
    elif [ "${REQUEST_STRING}" = "${WOODBURY_ACTION_EXECUTE}" ]; then
        echo "${WOODBURY_ACTION_EXECUTE}"
    else
        echo "Unknown action specified '${REQUEST_STRING}'!  Expected either '${WOODBURY_ACTION_PREPARE}' or '${WOODBURY_ACTION_EXECUTE}'." >&2
    fi
}
