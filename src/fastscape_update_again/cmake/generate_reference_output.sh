#!/bin/bash

# This script executes all tests and overwrites
# the reference output in the source directory for all failing
# outputs. Finally, a file changes.diff is generated in the build directory
# with the diffs. The current directory is assumed to be the build directory,
# while the path of this script needs to be inside the source cmake/
# directory.

echo "Overwriting test output with reference output ..."

SRC_PATH=`dirname $0`
SRC_PATH=`cd $SRC_PATH/..;pwd`
OUT=$PWD/changes.diff

if [ "$ASPECT_TESTS_VERBOSE" == "1" ]; then
  ASPECT_GENERATE_REFERENCE_OUTPUT=1 ctest -j 4
else
  ASPECT_GENERATE_REFERENCE_OUTPUT=1 ctest -j 4 >/dev/null
fi

cd $SRC_PATH
git diff tests/ >$OUT
if [ -s $OUT ]; then
  echo "generated patch file: $OUT";
  echo "modified files:"
  git diff --name-only tests/
else
  echo "no reference file changed."
fi
