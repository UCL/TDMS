# **`MATLAB` Benchmarking Scripts**

The directory `matlab_benchmark_scripts` contains scripts which perform band-limited interpolation (BLi) using `MATLAB`'s `interp` function.
`TDMS`'s interpolation schemes are based off this `MATLAB` function (specficially, in the coefficients the scheme uses to interpolate).

In order to test that the interpolation is correctly implimented in the `TDMS` source, we provide unit tests that benchmark against `MATLAB`'s implimentations. These scripts are provided here for developer use and documentation.

Currently; the unit tests have the output error values from `MATLAB` hard-coded into the source, and are required to achieve a comparable level of accuracy for the tests to pass. Each script contains a note referencing which unit test it provides a benchmark for.