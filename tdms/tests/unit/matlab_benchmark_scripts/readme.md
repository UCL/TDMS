# **`MATLAB` Benchmarking Scripts**

This directory contains scripts that generate MATLAB `.mat` files for use in the unit tests of TDMS, or provide benchmarks for the units that are tested.

### Unit test `.mat` data

A number of our unit tests require the presence of a `.mat` file to read/write data from/to during testing.
Running the following scripts in a MATLAB session within this directory will produce these `.mat` files.

### Band-limited Interpolation Benchmarking

The `benchmark_` scripts perform band-limited interpolation (BLi) using `MATLAB`'s `interp` function.
`TDMS`'s interpolation schemes are based off this `MATLAB` function (specficially, in the coefficients the scheme uses to interpolate).

In order to test that the interpolation is correctly implimented in the `TDMS` source, we provide unit tests that benchmark against `MATLAB`'s implimentations. These scripts are provided here for developer use and documentation.

Currently; the unit tests have the output error values from `MATLAB` hard-coded into the source, and are required to achieve a comparable level of accuracy for the tests to pass. Each script contains a note referencing which unit test it provides a benchmark for.
