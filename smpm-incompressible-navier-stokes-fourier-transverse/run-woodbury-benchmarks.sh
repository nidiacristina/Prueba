#!/bin/sh

# benchmark harness for the Woodbury Identity-based Navier-Stokes equation
# solver.  this is intended to explore the space of MPI ranks and OpenMP
# threads across one or more test cases and present a concise table of
# timing information to help guide design and optimization efforts.
#
# NOTE: this is only intended to be run on a single host.  since most
#       design related effects of communication may be seen within MPI
#       ranks on a single system, coupled with the fact that it is harder
#       to handle over subscription across multiple nodes, this script
#       does not attempt to handle benchmarks involving more than one
#       host.

# display the usage on standard output.
print_usage()
{
    echo "Usage: $0 [-h] [-L <launcher>] [-M <ranks>] [-O <threads>] [-T <tag>] <solver path> <viewer path> <test directory> [<test name> [...]]"
    echo
    echo "Benchmark harness for the direct solver for the Navier-Stokes equation using the Woodbury"
    echo "identity.  Executes the solver, at <solver path>, for one or more test cases found in"
    echo "<test directory>.  By default, all test cases found in <test directory> are executed,"
    echo "though specific tests may be run by providing <test name>'s on the command line.  Timing"
    echo "results are displayed on standard output, first grouped by test, then sorted according"
    echo "to performance, by using information produced by the SMPM Viewer found at <viewer path>."
    echo
    echo "Each test case is executed is invoked multiple times, iterating across multiple"
    echo "combinations of MPI ranks and OpenMP threads.  The space of ranks and thread counts"
    echo "may be restricted via the command line options -M and -O, or derived from the system's"
    echo "processor count and the input file of each test executed.  Additionally, the space"
    echo "benchmarked will be reduced, in either ranks or threads, or both, if it would exceed"
    echo "the system's total core count."
    echo
    echo "Each benchmarked test case's outputs will be named such that each combination of"
    echo "rank and thread count is unique.  The form of the name is:"
    echo
    echo "   <case name>_MPI=<rank count>_OMP=<thread count>_out[.stdout]"
    echo
    echo "The outputs may be optionally prefixed with a tag via the -T option."
    echo
    echo "The number of MPI ranks and OpenMP threads may be specified via the -M and -O options,"
    echo "respectively.  Both can specify a list of cores to use during benchmarking, either"
    echo "by explicitly enumerating, or by implicitly specifying ranges.  The argument accepted"
    echo "by both options takes the form of:"
    echo
    echo "     <range>[,<range>[...]]"
    echo
    echo "Where each <range> may either be a single number, or may be a range in the form N-M,"
    echo "where N is no greater than M."
    echo
    echo "The command line flags shown above are described below:"
    echo
    echo "   -h            Display this help message and exit."
    echo "   -L <launcher> Uses <launcher> to invoke the solver instead of ${MPI_LAUNCHER}."
    echo "   -M <ranks>    Specifies the number of MPI ranks to use during benchmarking, where"
    echo "                 <ranks> has the form of <range> above.  If omitted, defaults to all"
    echo "                 positive divisors greater than two for each test case benchmarked."
    echo "   -O <threads>  Specifies the OpenMP thread count to use during benchmarking, where"
    echo "                 <threads> has the form of <range> above.  If omitted, defaults to"
    echo "                 one to N, where N is the largest integer whose product with the"
    echo "                 smallest suitable MPI rank is smaller than the system's core count."
    echo "                 By default, all available cores will be used, tests permitting."
    echo "   -T <tag>      Use <tag> as a prefix for generated output files during test"
    echo "                 execution.  This allows successive executions to not overwrite"
    echo "                 each others output, and is particularly useful when test case"
    echo "                 debugging has been enabled.  By default, no tag is used."
    echo
}

# takes an executable, a path to the SMPM log viewer, and a test directory,
# and benchmarks each of the sub-tests within the specified directory.  it is
# assumed that the test case directory exists.  each sub-test is invoked
# multiple times, based on the execution environment described by MPI_RANKS
# and OPENMP_THREADS, restricted to the ranks and threads that can
# successfully execute each sub-test.  the timing information for each
# benchmark is summarized and displayed on standard output.
run_benchmark_case()
{
    EXECUTABLE=$1
    VIEWER=$2
    DIRECTORY=$3

    # move into the test directory.
    cd ${DIRECTORY}

    # identify all of the input files in this test directory.  there should
    # only be one, though this allows multiple to be aggregated within a
    # single directory.
    INPUT_FILES=`find ${DIRECTORY} -type f -name "*_in"`

    echo

    # walk through each input file and invoke the executable with it.
    for INPUT_FILE in ${INPUT_FILES}; do
        # identify the input and output files necessary for execution and
        # verification.
        RELATIVE_INPUT_FILE=`basename ${INPUT_FILE}`
        RELATIVE_OUTPUT_BASE=`echo ${RELATIVE_INPUT_FILE} | sed -e 's#_in$#_out#'`

        # if we have a tag, prefix the generated output file with it.
        if [ -n "${OUTPUT_TAG}" ]; then
            RELATIVE_OUTPUT_BASE="${OUTPUT_TAG}-${RELATIVE_OUTPUT_BASE}"
        fi

        # create a short name for the test that can be used for communicating
        # success/failure when we're appropriately verbose.
        TEST_NAME=`echo ${RELATIVE_INPUT_FILE} | sed -e 's#_in$##'`

        # announce this input file within this test case.
        echo "  ${TEST_NAME}"
        echo
        echo "${BENCHMARK_HEADER_FORMAT}"

        # reset the benchmark output for this input file.
        BENCHMARK_OUTPUT=

        # generate the MPI ranks we're interested in for this test case.
        CANDIDATE_MPI_RANKS=`get_valid_mpi_ranks ${RELATIVE_INPUT_FILE} ${MPI_RANKS}`

        # loop over the MPI ranks of interest.
        for MPI_RANK in ${CANDIDATE_MPI_RANKS}; do
            # generate the OpenMP thread counts that make sense for the
            # current MPI rank.
            CANDIDATE_OPENMP_THREADS=`get_valid_openmp_threads ${MPI_RANK} ${OPENMP_THREADS}`

            # iterate over the OpenMP thread counts of interest.
            for OPENMP_THREAD in ${CANDIDATE_OPENMP_THREADS}; do
                # add the current MPI rank and OpenMP thread count into the
                # output file name so that we can sort things out after the
                # fact.
                RELATIVE_OUTPUT_FILE=`echo ${RELATIVE_OUTPUT_BASE} | sed -e s/_out/_MPI=${MPI_RANK},OMP=${OPENMP_THREAD}_out/`
                RELATIVE_CONSOLE_OUTPUT_FILE=${RELATIVE_OUTPUT_FILE}.stdout

                # remove the previous run's log file.  if the solver fails to execute
                # successfully or does not modify the log, our execution results are
                # incorrect.
                rm -f ${RELATIVE_OUTPUT_FILE}

                # run the solver and store the output.
                OMP_NUM_THREADS=${OPENMP_THREAD} ${MPI_LAUNCHER} -n ${MPI_RANK} ${EXECUTABLE} ${RELATIVE_INPUT_FILE} ${RELATIVE_OUTPUT_FILE} > ${RELATIVE_CONSOLE_OUTPUT_FILE}
                EXIT_STATUS=$?

                # see if something went awry and let the user know.
                if [ ! -f ${RELATIVE_OUTPUT_FILE} ]; then
                    echo "Simulation did not create ${RELATIVE_OUTPUT_FILE}." >&2
                    continue
                elif [ ${EXIT_STATUS} -ne 0 ]; then
                    echo "Simulation returned with a non-zero status of ${EXIT_STATUS} (MPI=${MPI_RANK}, OMP=${OPENMP_THREAD})." >&2
                    continue
                fi

                # pull out the results from this run.
                #
                # NOTE: this assumes that runs take more than one timestep so
                #       that the total wall-time is greater than the timestep
                #       time.
                #
                # NOTE: eventually we should parse the output of the log file
                #       via the viewer at ${VIEWER_PATH}.
                TIMINGS=`egrep -e 'Wall time (elapsed|per time-step)' ${RELATIVE_CONSOLE_OUTPUT_FILE} | \
                         awk '{ print $(NF-1) }' | sort -n | xargs`

                TIME_TIMESTEP=`echo ${TIMINGS} | awk '{ print $1 }'`
                TIME_TOTAL=`echo ${TIMINGS} | awk '{ print $2 }'`

                BENCHMARK_OUTPUT="${BENCHMARK_OUTPUT} ${MPI_RANK}:${OPENMP_THREAD}:${TIME_TIMESTEP}:${TIME_TOTAL}"
            done
        done

        # produce formatted benchmark results and then sort by the desired
        # column.
        format_benchmark_results "${BENCHMARK_OUTPUT}" | \
            sort -k ${BENCHMARK_RESULT_SORT_FIELD}

        echo
    done
}

# takes a space-delimited list of colon-delimited benchmark results and pretty
# prints them to standard output.
format_benchmark_results()
{
    BENCHMARK_RESULTS="$1"

    # walk through each of the benchmarks and pretty prinnt its results.
    for RESULT in ${BENCHMARK_RESULTS}; do
        MPI_RANK_COUNT=`echo ${RESULT} | cut -f1 -d:`
        OPENMP_THREAD_COUNT=`echo ${RESULT} | cut -f2 -d:`
        TIME_TIMESTEP=`echo ${RESULT} | cut -f3 -d:`
        TIME_TOTAL=`echo ${RESULT} | cut -f4 -d:`

        # compute the total core count used for this benchmark and use it
        # to compute a normalized time.
        BENCHMARK_CORE_COUNT=`expr ${MPI_RANK_COUNT} '*' ${OPENMP_THREAD_COUNT}`
        TIME_PER_CORE=`echo ${TIME_TOTAL} / ${BENCHMARK_CORE_COUNT} | bc -l`

        printf "${BENCHMARK_RESULT_FORMAT}\n" \
            ${MPI_RANK_COUNT} \
            ${OPENMP_THREAD_COUNT} \
            ${TIME_TIMESTEP} \
            ${TIME_TOTAL} \
            ${TIME_PER_CORE}
    done
}

# takes a solver input file name and a list of space delimited, potential MPI
# ranks to run the solver with, and echoes the ranks which could execute the
# solver without error.  order within the supplied potential ranks is
# preserved within the ranks output.
get_valid_mpi_ranks()
{
    INPUT_FILENAME=$1
    shift 1

    # everything else provided is a potential rank to evaluate.
    POTENTIAL_MPI_RANKS="$*"

    # start out with an empty list of ranks that are valid for the input file.
    VALID_MPI_RANKS=

    # the solver can only distribute work across MPI ranks in x-subdomain
    # sized quantitized, so determine how many of those we have.
    NUMBER_SUBDOMAINS=`get_number_x_subdomains ${INPUT_FILENAME}`

    # walk through each potential rank and compute a suitable core count from
    # it.  if any other value is returned, then that particular number of MPI
    # ranks could not evenly distribute the workload and is not a valid
    # candidate.
    for POTENTIAL_MPI_RANK in ${POTENTIAL_MPI_RANKS}; do
        COMPUTED_MPI_RANK=`compute_suitable_core_count ${POTENTIAL_MPI_RANK} ${NUMBER_SUBDOMAINS}`

        if [ ${COMPUTED_MPI_RANK} -eq ${POTENTIAL_MPI_RANK} ]; then
            VALID_MPI_RANKS="${VALID_MPI_RANKS} ${POTENTIAL_MPI_RANK}"
        fi
    done

    echo ${VALID_MPI_RANKS}
}

# takes the number of MPI ranks and a list of space delimited, potential
# OpenMP thread counts to run the solver with, and echoes the thread counts
# which could execute them solver without oversubscribing the system.  order
# within the supplied potential threads is preserved within the output.
get_valid_openmp_threads()
{
    REFERENCE_MPI_RANK=$1
    shift

    # everything else provided is a potential thread count to evaluate.
    POTENTIAL_OPENMP_THREADS="$*"

    # start out with an empty list of thread counts that are valid with
    # respect to the reference MPI rank count.
    VALID_OPENMP_THREADS=

    # walk through each potential thread count and ensure this does not
    # exceed the total system core count.  if it doesn't, then it is a
    # valid thread count.
    for POTENTIAL_OPENMP_THREAD in ${POTENTIAL_OPENMP_THREADS}; do
        if [ `expr ${REFERENCE_MPI_RANK} '*' ${POTENTIAL_OPENMP_THREAD}` -le ${SYSTEM_CORE_COUNT} ]; then
            VALID_OPENMP_THREADS="${VALID_OPENMP_THREADS} ${POTENTIAL_OPENMP_THREAD}"
        fi
    done

    echo ${VALID_OPENMP_THREADS}
}

# takes a comma-delimited list of positive numbers, or ranges of positive
# numbers, and outputs either "yes" or "no" depending on whether the range is
# considered valid or not, respectively.  ranges are a separated by a dash and
# must be in non-decreasing order.
is_valid_range_set()
{
    INPUT_RANGES="$1"

    # valid input must be one of the following:
    #
    #  * a positive integer.
    #  * a range of positive integers, separated by a dash (-), in
    #    non-decreasing order.
    #
    # convert commas into spaces and validate the individual tokens.
    #
    # NOTE: this does not ensure that the input is precisely the above, but
    #       is (likely) good enough for what we'll encounter.
    for INPUT_TOKEN in `echo ${INPUT_RANGES} | tr , ' '`; do
        echo ${INPUT_TOKEN} | egrep -q -e '^([0-9]+|[0-9]+-[0-9]+)$'
        if [ $? -ne 0 ]; then
            echo "no"
            return
        fi

        # if this is a range, ensure that it is not decreasing.
        echo ${INPUT_TOKEN} | grep -q '-'
        if [ $? -eq 0 ]; then
            RANGE_START=`echo ${INPUT_TOKEN} | awk -F - '{ print $1 }'`
            RANGE_STOP=`echo ${INPUT_TOKEN} | awk -F - '{ print $2 }'`

            if [ ${RANGE_START} -gt ${RANGE_STOP} ]; then
                echo "no"
                return
            fi
        fi

        # this token looks sensible.  move to the next one.
    done

    echo "yes"
}

# takes a comma-delimited list of positive numbers, or ranges of positive
# numbers, and outputs a space-delimited sortued, unique list of numbers
# that were contained in the input.
expand_ranges()
{
    # NOTE: we assume our input has been validated.  see is_valid_range_set()
    #       for a description of what valid means.
    INPUT_RANGES="$1"

    OUTPUT_LIST=

    for INPUT_TOKEN in `echo ${INPUT_RANGES} | tr , ' '`; do
        # if this is a range, add in all numbers in the range.  otherwise,
        # add the token itself.
        echo ${INPUT_TOKEN} | grep -q '-'
        if [ $? -eq 0 ]; then
            RANGE_START=`echo ${INPUT_TOKEN} | awk -F - '{ print $1 }'`
            RANGE_STOP=`echo ${INPUT_TOKEN} | awk -F - '{ print $2 }'`

            OUTPUT_LIST="${OUTPUT_LIST} `seq ${RANGE_START} ${RANGE_STOP}`"
        else
            OUTPUT_LIST="${OUTPUT_LIST} ${INPUT_TOKEN}"
        fi
    done

    # output a sorted list with duplicates removed.
    echo ${OUTPUT_LIST} | sort -u | xargs
}

# pull in our utility functions.
. `dirname $0`/woodbury-utility.sh

# prefix tag for generated files.  this allows for multiple invocations of the
# test framework to be run, each with a different tag, without collisions
# between the files generated.  this is typically provided in conjunction with
# debug output to prevent files from being deleted after execution.
OUTPUT_TAG=

# launcher used for MPI applications.  most MPI distributions use the default,
# mpirun, though certain sites (e.g. TACC) have different launchers that are
# aware of each system's topology.
MPI_LAUNCHER=mpirun

# compute the total number of cores in the system.  this is used to compute
# the threshold for over-subscription as well as the upper bounds for the
# maximum number of MPI ranks and OpenMP threads.
SYSTEM_CORE_COUNT=`get_system_core_count`

# header and result format string for pretty printing benchmark results.
# these are sized with the assumption that test cases will not run for more
# than 2-3 hours, and will live on a single node.
BENCHMARK_HEADER_FORMAT="        MPI   OpenMP  sec/step   Total     Time/Core"
BENCHMARK_RESULT_FORMAT="        %2d    %3d    %8.5f    %7.4f   %8.5f"

# sort benchmark results within a single test by the time per timestep field.
BENCHMARK_RESULT_SORT_FIELD=3

# parse the command line options provided.
while getopts "hLM:O:T:" FLAGS
do
  case ${FLAGS} in
      h)
          print_usage
          exit 0
          ;;
      L)
          MPI_LAUNCHER=${OPTARG}
          ;;
      M)
          MPI_RANKS=${OPTARG}
          ;;
      O)
          OPENMP_THREADS=${OPTARG}
          ;;
      T)
          OUTPUT_TAG=${OPTARG}
          ;;
  esac
done

# skip over any command line options provided.
shift `expr ${OPTIND} - 1`

# ensure that we got at least the minimum number of arguments necessary to
# operate.
if [ $# -lt 3 ]; then
    echo "Incorrect number of arguments."
    print_usage

    exit 1
fi

# identify the solver, the viewer, and the location of its test cases.
SOLVER_PATH=$1
VIEWER_PATH=$2
TESTING_DIRECTORY=$3

# skip past the paths and the test directory.  everything remaining is a test
# name that the user has requested to run.
shift 3
TEST_CASE_CANDIDATES="$@"

# verify that the solver exists and the testing directory is present.
if [ ! -x "${SOLVER_PATH}" ]; then
    echo "The supplied solver (${SOLVER_PATH}) is not executable!" >&2

    exit 1
elif [ ! -x "${VIEWER_PATH}" ]; then
    echo "The supplied viewer (${VIEWER_PATH}) is not executable!" >&2

    exit 1
elif [ ! -d "${TESTING_DIRECTORY}" ]; then
    echo "The testing directory (${TESTING_DIRECTORY}) is not a directory!" >&2

    exit 1
fi

# ensure that the MPI launcher is executable.
MPI_LAUNCHER=`which "${MPI_LAUNCHER}" 2>/dev/null`
if [ ! -x "${MPI_LAUNCHER}" ]; then
    echo "${MPI_LAUNCHER} is not executable!" >&2
fi

# ensure that the ranks requested make sense, and can be expanded into
# something usable.
if [ -n "${MPI_RANKS}" ]; then
    if [ `is_valid_range_set ${MPI_RANKS}` != "yes" ]; then
        echo "Invalid MPI ranks specified (${MPI_RANKS})." >&2

        exit 1
    fi

    # expand any compressed ranges into a space delimited list.
    MPI_RANKS=`expand_ranges ${MPI_RANKS}`
else
    # the user doesn't have specific preferences on which ranks to benchmark.
    # default to all possible ranks.
    MPI_RANKS=`seq 2 ${SYSTEM_CORE_COUNT}`
fi

# ensure that the threads requested make sense, and can be expanded into
# something usable.
if [ -n "${OPENMP_THREADS}" ]; then
    if [ `is_valid_range_set ${OPENMP_THREADS}` != "yes" ]; then
        echo "Invalid OpenMP thread count specified (${OPENMP_THREADS})." >&2

        exit 1
    fi

    # expand any compressed ranges into a space delimited list.
    OPENMP_THREADS=`expand_ranges ${OPENMP_THREADS}`
else
    # the user doesn't have specific preference on how many threads to
    # benchmark.  default to a range covering serial to fully parallel
    # on this system.
    #
    # NOTE: this over estimates the thread counts that are viable since
    #       a single MPI rank is not currently allowed (2014/06/25).
    #       rather than try to correctly identify the upper bound, we'll
    #       let it be pruned down within each test case.
    OPENMP_THREADS=`seq 1 ${SYSTEM_CORE_COUNT}`
fi

# normalize the solver paths and top level test directory so that we do not have
# to be careful of which directory we're in while running tests.
SOLVER_PATH=`readlink -f ${SOLVER_PATH}`
VIEWER_PATH=`readlink -f ${VIEWER_PATH}`
TESTING_DIRECTORY=`readlink -f ${TESTING_DIRECTORY}`

# find all of the available test cases beneath the testing directory.
TEST_CASES_AVAILABLE=`enumerate_available_test_cases "${TESTING_DIRECTORY}"`

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

    run_benchmark_case "${SOLVER_PATH}" "${VIEWER_PATH}" "${TESTING_DIRECTORY}/${TEST_CASE}"

    increment_variable TEST_NUMBER
done
