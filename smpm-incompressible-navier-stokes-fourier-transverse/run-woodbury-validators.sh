#!/bin/sh

# test harness for the SMPM Validator.

# display the usage on standard output.
print_usage()
{
    echo "Usage: $0 [-A <actions>] [-h] [-L <launcher>] [-M <path>] [-O] [-p <cores>] <validator path> <validation directory> [<test name> [...]]"
    echo
    echo "Test harness for the SMPM Validator.  Executes the validator at <solve path>,"
    echo "for one or more test cases found in <validation directory>.  By default,"
    echo "all test cases found in <validation directory> are executed, though specific test may"
    echo "be run by providing <test name>'s on the command line."
    echo
    echo "Each test case's metadata specifies an accuracy threshold and a list of validation"
    echo "routines for the SMPM Validator to execute.  Test cases pass if each routine meets"
    echo "or exceeds the required threshold, and the results of each test are printed to"
    echo "standard output.  The test harness returns success if all test cases requested"
    echo "pass, otherwise returns failure."
    echo
    echo "The command line flags shown above are described below:"
    echo
    echo "   -A <actions>  A list of actions specifying the validator should be run in stages"
    echo "                 instead of in their entirety.  This allows preparation to be performed"
    echo "                 on one machine and execution on another, so as to minimize the"
    echo "                 dependencies on the execution environment.  <actions> must be"
    echo "                 one of \"${WOODBURY_ACTION_PREPARE}\" or \"${WOODBURY_ACTION_EXECUTE}\".  If <actions> is omitted it defaults"
    echo "                 to both of them."
    echo "   -h            Display this help message and exit."
    echo "   -L <launcher> Uses <launcher> to invoke the solver instead of ${MPI_LAUNCHER}."
    echo "   -M <path>     Specifies <path> as the location of the SMPM MATLAB code.  If omitted"
    echo "                 defaults to ${MATLAB_PATH}."
    echo "   -O            Indicates that Octave should be used instead of MATLAB for test"
    echo "                 preparation.  Otherwise Octave will only be used if MATLAB cannot be"
    echo "                 found on the system."
    echo "   -p <cores>    Specifies the maximum number of cores to use while executing test"
    echo "                 cases.  Note that fewer than <cores> many processors may be used"
    echo "                 for individual test cases, should the test's grid not support it."
    echo "                 By default, all available cores will be used, test configuration"
    echo "                 permitting."
    echo
}

# takes a validator executable, a validation test directory, and an input file
# and runs the validator on each of the tests found in the test directory.  if
# each of the tests passes with the accuracy specified in the test directory,
# success is recorded, failure otherwise.  results are recorded for later
# reporting via report_results().
run_test_case()
{
    VALIDATOR=$1
    DIRECTORY=$2
    INPUT_FILE=$3

    # we need to move into the test directory since the paths we setup below
    # are all relative.
    cd ${DIRECTORY}

    # identify the input and output files necessary for execution and
    # verification.
    RELATIVE_INPUT_PREFIX="`basename ${INPUT_FILE} | sed -e 's#_in##'`"

    RELATIVE_INPUT_FILE="${RELATIVE_INPUT_PREFIX}_in"
    RELATIVE_CONSOLE_OUTPUT_FILE="${RELATIVE_INPUT_PREFIX}_log"

    RELATIVE_METADATA_FILE="validation.txt"

    # create a short name for the test that can be used for communicating
    # success/failure when we're appropriately verbose.
    TEST_NAME="${RELATIVE_INPUT_PREFIX}"

    # if we do not have metadata, we can't run the validator.
    if [ ! -f ${RELATIVE_METADATA_FILE} ]; then
        test_failed ${TEST_NAME} "${RELATIVE_METADATA_FILE} is missing"
        return
    fi

    # get the parameters for the validator.
    TEST_TOLERANCE=`head -n1 ${RELATIVE_METADATA_FILE}`
    TEST_ROUTINES=`sed -n -e '2,$p' ${RELATIVE_METADATA_FILE}`

    # determine how many of the cores we've been allocated should be used
    # to execute this test.
    #
    # XXX: this needs to go into the input generator to make sure the
    #      preparation and execution phases are synchronized.
    #
    NUMBER_PROCESSORS_TEST=`get_suitable_core_count ${NUMBER_PROCESSORS} ${RELATIVE_INPUT_FILE}`
    NUMBER_THREADS_TEST=`get_largest_thread_count ${NUMBER_PROCESSORS_TEST} ${NUMBER_PROCESSORS}`

    # run the validator and store the output.
    OMP_NUM_THREADS=${NUMBER_THREADS_TEST} ${MPI_LAUNCHER} -n ${NUMBER_PROCESSORS_TEST} ${VALIDATOR} ${RELATIVE_INPUT_FILE} ${TEST_TOLERANCE} ${TEST_ROUTINES} >${RELATIVE_CONSOLE_OUTPUT_FILE} 2>&1
    EXIT_STATUS=$?

    if [ ${EXIT_STATUS} -eq 0 ] ; then
        test_passed ${TEST_NAME}
    else
        test_failed ${TEST_NAME} "${FAILURE_MESSAGE}"
    fi
}

# takes a number of successes and failures and displays the them, along with
# the relative percentages, to standard output.
report_results()
{
    NUMBER_SUCCESSES=${1:-0}
    NUMBER_FAILURES=${2:-0}

    NUMBER_TESTS=`expr ${NUMBER_SUCCESSES} + ${NUMBER_FAILURES}`

    echo
    echo "Tests executed: ${NUMBER_TESTS}"
    echo

    if [ ${NUMBER_TESTS} != 0 ]; then
        PERCENTAGE_SUCCESSES=`echo ${NUMBER_SUCCESSES} / ${NUMBER_TESTS} '*' 100.0 | bc -l`
        PERCENTAGE_FAILURES=`echo ${NUMBER_FAILURES} / ${NUMBER_TESTS} '*' 100.0 | bc -l`
        printf "     Successes: %3d    (%6.2f %%)\n" ${NUMBER_SUCCESSES} ${PERCENTAGE_SUCCESSES}
        printf "      Failures: %3d    (%6.2f %%)\n" ${NUMBER_FAILURES} ${PERCENTAGE_FAILURES}
        echo
    fi
}

# takes the number of MPI ranks and a list of space delimited, potential
# OpenMP thread counts to run the solver with, and echoes the thread counts
# which could execute them solver without oversubscribing the system.  order
# within the supplied potential threads is preserved within the output.
get_largest_thread_count()
{
    REFERENCE_MPI_RANK=$1
    OVERSUBSCRIPTION_THRESHOLD=$2

    THREAD_COUNT=1

    # walk through each potential thread count and ensure this does not
    # exceed the total system core count.  if it doesn't, then it is a
    # valid thread count.
    for THREAD_COUNT in `seq ${OVERSUBSCRIPTION_THRESHOLD} -1 1`; do
        TOTAL_THREAD_COUNT=`expr ${REFERENCE_MPI_RANK} '*' ${THREAD_COUNT}`

        if [ `expr ${REFERENCE_MPI_RANK} '*' ${THREAD_COUNT}` -le ${OVERSUBSCRIPTION_THRESHOLD} ]; then
            break
        fi
    done

    echo ${THREAD_COUNT}
}

# records the specified named test as a success.  a message, possibly indented
# using a caller supplied string, is displayed to standard output indicating
# as such.  the global variable SUCCESS_COUNT is incremented by one.
test_passed()
{
    NAME=$1
    MESSAGE=${2:-""}
    INDENTATION=${3:-"    "}

    increment_variable SUCCESS_COUNT

    if [ -n "${MESSAGE}" ]; then
        MESSAGE=" (${MESSAGE})"
    fi

    echo "${INDENTATION}Test \"${NAME}\" - Passed${MESSAGE}"
}

# records the specified named test as a failure.  a message, possibly indented
# using a caller supplied string, is displayed to standard output indicating
# as such.  the global variable FAILURE_COUNT is incremented by one.
test_failed()
{
    NAME=$1
    MESSAGE=${2:-""}
    INDENTATION=${3:-"    "}

    increment_variable FAILURE_COUNT

    if [ -n "${MESSAGE}" ]; then
        MESSAGE=" (${MESSAGE})"
    fi

    echo "${INDENTATION}Test \"${NAME}\" - Failed${MESSAGE}"
}

# pull in our utility functions.
. `dirname $0`/woodbury-utility.sh

# counters indicate how many tests have passed and failed.  these are
# manipulated through the test_passed()/test_failed() functions as individual
# tests are executed.
SUCCESS_COUNT=0
FAILURE_COUNT=0

# launcher used for MPI applications.  most MPI distributions use the default,
# mpirun, though certain sites (e.g. TACC) have different launchers that are
# aware of each system's topology.
MPI_LAUNCHER=mpirun

# path to the SMPM MATLAB utilities.
MATLAB_PATH="./smpm_matlab_utilities"

# our test case input generators are MATLAB functions.  we attempt to use MATLAB
# to invoke them unless told not to try.  specifying this on the command line
# simply skips the detection of a working MATLAB installation.
OCTAVE_FLAG=no

# the number of processors that the solver is invoked with.  by default, use
# all of the available resources on the system.
NUMBER_PROCESSORS=`get_system_core_count`

# the actions the validation framework should perform by default.  specifying a
# subset of available actions, such as to only prepare or execute previously
# prepared cases, can be specified on the command line.
HARNESS_ACTIONS=`parse_test_actions`

# parse the command line options provided.
while getopts "A:hL:M:Op:" FLAGS
do
  case ${FLAGS} in
      A)
          # get a list of actions.  we assume the parser emits error messages
          # when invalid input is supplied.
          HARNESS_ACTIONS=`parse_test_actions ${OPTARG}`
          if [ -z "${HARNESS_ACTIONS}" ]; then
              exit 1
          fi
          ;;
      h)
          print_usage
          exit 0
          ;;
      L)
          MPI_LAUNCHER=${OPTARG}
          ;;
      M)
          MATLAB_PATH=${OPTARG}
          ;;
      O)
          OCTAVE_FLAG=yes
          ;;
      p)
          NUMBER_PROCESSORS=${OPTARG}
          ;;
  esac
done

# skip over any command line options provided.
shift `expr ${OPTIND} - 1`

# ensure that we got at least the minimum number of arguments necessary to
# operate.
if [ $# -lt 2 ]; then
    echo "Incorrect number of arguments."
    print_usage

    exit 1
fi

# identify the validator and the location of its test cases.
VALIDATOR_PATH=$1
VALIDATION_DIRECTORY=$2

# map our list of actions into flags to enable/disable execution paths.
PREPARE_FLAG=`echo ${HARNESS_ACTIONS} | grep -q ${WOODBURY_ACTION_PREPARE} && echo yes || echo no`
EXECUTE_FLAG=`echo ${HARNESS_ACTIONS} | grep -q ${WOODBURY_ACTION_EXECUTE} && echo yes || echo no`

# skip past the path and the test directory.  everything remaining is a test
# name that the user has requested to run.
shift 2
TEST_CASE_CANDIDATES="$@"

# verify that the validator exists, and that the testing directory is present.
if [ ! -x "${VALIDATOR_PATH}" ]; then
    echo "The supplied validator (${VALIDATOR_PATH}) is not executable!" >&2

    exit 1
elif [ ! -d "${VALIDATION_DIRECTORY}" ]; then
    echo "The validation directory (${VALIDATION_DIRECTORY}) is not a directory!" >&2

    exit 1
fi

# find something that can interpret our MATLAB code.  take care to respect
# the caller's wishes for using Octave.
MATLAB_COMMAND=`get_matlab_command ${OCTAVE_FLAG}`

if [ -z "${MATLAB_COMMAND}" ]; then
    echo "Could not find something to run MATLAB code!" >&2

    exit 1
fi

# ensure that the number of processors requested is valid.  since this is
# independent of any specific tests we'll run, this essentially checks that it
# is a positive number.
if [ `is_valid_processor_count ${NUMBER_PROCESSORS}` != "yes" ]; then
    echo "Invalid processor count (${NUMBER_PROCESSORS}) specified!" >&2

    exit 1
fi

# ensure that the MPI launcher is executable.
MPI_LAUNCHER=`which "${MPI_LAUNCHER}" 2>/dev/null`
if [ ! -x "${MPI_LAUNCHER}" ]; then
    echo "${MPI_LAUNCHER} is not executable!" >&2
fi

# normalize the important paths so that we do not have to be careful of which
# directory we're in while running tests.
VALIDATOR_PATH=`readlink -f ${VALIDATOR_PATH}`
MATLAB_PATH=`readlink -f ${MATLAB_PATH}`
VALIDATION_DIRECTORY=`readlink -f ${VALIDATION_DIRECTORY}`

# find all of the available test cases beneath the testing directory.
TEST_CASES_AVAILABLE=`enumerate_available_test_cases "${VALIDATION_DIRECTORY}"`

# set the list of test cases to run.  if the user requested specific tests to
# run, use those.  otherwise, use everything available.
if [ -z "${TEST_CASE_CANDIDATES}" ]; then
    TEST_CASE_CANDIDATES="${TEST_CASES_AVAILABLE}"
fi

TEST_CASES=`select_test_cases "${TEST_CASES_AVAILABLE}" "${TEST_CASE_CANDIDATES}"`

# ensure that we have something to do before proceeding.
if [ -z "${TEST_CASES}" ]; then
    echo "None of the requested test cases are valid.  Exiting." >&2

    exit 1
fi

# keep track of which test we're on.
TEST_NUMBER=1

# iterate through each test and run it.
for TEST_CASE in ${TEST_CASES}; do
    # add a little whitespace between successive test case outputs after the
    # first.
    if [ ${TEST_NUMBER} -gt 1 ]; then
        echo
    fi

    # announce the test case and its position relative to the start of these
    # tests.
    echo "#${TEST_NUMBER} - Test case ${TEST_CASE}: "

    # setup and execute this test case.
    #
    # NOTE: since we need this function to modify global state (read: test
    #       failures) we cannot wrap this in backticks and evaluate it in a
    #       subshell.  soooo, we have to result to some *nasty* coupling where
    #       prepare_test_case directly sets INPUT_FILE so we can use it
    #       below.
    #
    prepare_test_case "${VALIDATION_DIRECTORY}/${TEST_CASE}" ${PREPARE_FLAG}

    if [ -n "${INPUT_FILE}" ]; then
        if [ ${EXECUTE_FLAG} = "yes" ]; then
            run_test_case "${VALIDATOR_PATH}" "${VALIDATION_DIRECTORY}/${TEST_CASE}" "${INPUT_FILE}" "${OUTPUT_TAG}"
        else
            echo "Skipped"
        fi
    fi

    increment_variable TEST_NUMBER
done

# report the results.
report_results ${SUCCESS_COUNT} ${FAILURE_COUNT}

# allow higher level logic to know how we did.
exit ${FAILURE_COUNT}
