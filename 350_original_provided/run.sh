# Due to MATLAB namespacing (or lack of concept thereof)
# this script only runs when the other files within this
# folder are present in the same directory as this script.

# Fail-fast on non-zero exit code
set -e

# Identify TDMS executable
if [ $# -gt 0 ]; then
    TDMS_EXE=$1
else
    TDMS_EXE="../build/tdms"
fi
echo "Looking for tdms executable at ${TDMS_EXE}"
TDMS_VERSION=$(${TDMS_EXE} --version)
echo "Found TDMS version ${TDMS_VERSION}"
echo "---"

# Name of the input file to generate, then pass to TDMS
if [ $# -gt 1 ]; then
    INPUT_FILE=$2
else
    INPUT_FILE="oct_ml_in"
fi

# Attempt to generate the input file
echo "Launching MATLAB..."

matlab -batch generate_input_file input_file_ml.m $INPUT_FILE

echo "---"
if [ $? -ne 0 ]; then
    echo "MATLAB could not generate the input file (${INPUT_FILE})"
    exit 1
fi

# Attempt to run TDMS
echo "Launching TDMS..."

${TDMS_EXE} ${INPUT_FILE}.mat 350_output.mat
