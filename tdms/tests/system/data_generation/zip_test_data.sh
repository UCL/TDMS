#!/bin/bash

# Check that we are in the data_generation folder
CURRENT_FOLDER=${PWD##*/}
CURRENT_FOLDER=${CURRENT_FOLDER:-/}

if [ ${CURRENT_FOLDER} != "data_generation" ]
then
    echo "Error: not in data_generation folder"
    exit 1
fi

if ! [ -d "${PWD}/../data" ]
then
    echo "Error: data folder does not exist - cannot place generated outputs there!"
    exit 1
fi

# zip_test_folder TEST_ID
function zip_test_folder
{
    # interpret inputs into actually readable things
    TEST_ID=$1
    TEST_DIR="arc_${TEST_ID}"

    # change into test directory
    OLD_WD=$PWD
    cd ${TEST_DIR}

    # create a zip folder in the data_generation directory for this test
    echo "=== Creating ${TEST_DIR}.zip"
    zip -rj "../../data/${TEST_DIR}.zip" .

    # return to top-level directory
    cd $OLD_WD
}

# start regenerating the test data
echo "Zipping TDMS system test data..."

zip_test_folder 01
zip_test_folder 02
zip_test_folder 03
zip_test_folder 08
zip_test_folder 09
zip_test_folder 10
zip_test_folder 12
zip_test_folder 13
zip_test_folder example_fdtd

echo "Done"
