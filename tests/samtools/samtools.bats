#!/bin/bash

load '/mnt/c/Users/abc/Desktop/WGS/tests/bats-support/load.bash'
load '/mnt/c/Users/abc/Desktop/WGS/tests/bats-assert/load.bash'
load '/mnt/c/Users/abc/Desktop/WGS/tests/bats-file/load.bash'

# create some variables for tests to use later
setup_file() {
  # a timestamp directory to contain test results
  datetime=$(date --iso-8601=seconds)
  test_dir=$PWD
  test_time_dir=$test_dir/$datetime
  mkdir $test_time_dir

#  input_bam = $test_dir/Inputs/example.bam
#  expected_outdir = $test_dir/Outputs

  # export these variables for use in other tests
  export test_time_dir
  export test_dir
#  export expected_outdir
#  export input_bam
}

setup() {
  mkdir $test_time_dir/${BATS_TEST_NUMBER}_$BATS_TEST_DESCRIPTION
  mkdir $test_time_dir/${BATS_TEST_NUMBER}_$BATS_TEST_DESCRIPTION/results
  outdir=$test_time_dir/${BATS_TEST_NUMBER}_$BATS_TEST_DESCRIPTION/results/
  export outdir
}

teardown() {
  sleep 1
}

# Run nextflow version
@test "RunNextflowVersion" {
  run nextflow -version # assert exit status is 0
  [ "${status}" -eq 0 ]
  echo "${output}" >$outdir/RunNextflowVersion.txt

  assert_file_not_empty $outdir/RunNextflowVersion.txt
}
