#!/bin/bash

# Check that MATLAB exists on the path
if ! [ -x "$(command -v matlab)" ]; then
  echo "Error: MATLAB is not on the path." >&2
  exit 1
fi

# Check that we are in the data_generation folder
CURRENT_FOLDER=${PWD##*/}
CURRENT_FOLDER=${CURRENT_FOLDER:-/}

if [ ${CURRENT_FOLDER} != "data_generation" ]
then
    echo "Error: not in data_generation folder"
    exit 1
fi

# regenerate_test_input TEST_ID BSCAN_FILE
function regenerate_test_input
{
    # interpret inputs into actually readable things
    TEST_ID=$1
    TEST_DIR="arc_${TEST_ID}"
    BSCAN_FILE=$2

    # change into test directory
    OLD_WD=$PWD
    cd ${TEST_DIR}

    # run bscan script
    matlab -nodisplay -r "${BSCAN_FILE} ; exit"

    # return to top-level directory
    cd $OLD_WD
}

regenerate_test_input 01 "run_pstd_bscan"
regenerate_test_input 02 "run_pstd_bscan"
regenerate_test_input 03 "run_pstd_bscan"
regenerate_test_input 08 "run_pstd_bscan_fast"
#regenerate_test_input 09 "run_pstd_bscan"
#regenerate_test_input 10 "run_pstd_bscan"
regenerate_test_input 12 "run_fdtd_bscan"
regenerate_test_input 13 "run_fdtd_bscan"
regenerate_test_input example_fdtd "run_fdtd_bscan"
