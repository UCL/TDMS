### MATLAB tests

The unit tests in this directory can be run with

```bash
matlab -batch run_tests
```

## `test_iteratefdtd_matrix_source_combinations`

The `iteratefdtd_matrix` function has a number of ways to produce the source terms from the information provided to it in an `input_file.m`.
Some combinations of the variables that _can_ be provided however are incompatible - either not providing sufficient information to generate the source fields, or providing contradictory information.

As such, the unit test `test_iteratefdtd_matrix_source_combinations` (provided as `arc_19` in [#208](https://github.com/UCL/TDMS/issues/208)) explicitly works through each possible combination of source-related inputs and checks for the corresponding error (or checks that no error is raised in the event the combination of inputs is valid).
The various combinations and expected results are listed in the table below.

`TD-field` is also know as `exi/eyi` in some docstrings.

| TD-field | usecd | compactsource | efname    | hfname    | Raises error? | Error info |
|:--------:|:-----:|:-------------:|:---------:|:---------:|:-------------:|:----------:|
| 1        | 1     | 1             | 1         | 1         | 1             | Cannot specify hfname if compact source |
| 1        | 1     | 1             | 1         | 0         | 0             |  |
| 1        | 1     | 1             | 0         | 1         | 1             | Cannot specify hfname if compact source |
| 1        | 1     | 1             | 0         | 0         | 0             |  |
| 1        | 1     | 0             | 1         | 1         | 0             |  |
| 1        | 1     | 0             | 1         | 0         | 1             | If not compact source, both efname and hfname must be specified |
| 1        | 1     | 0             | 0         | 1         | 1             | If not compact source, both efname and hfname must be specified |
| 1        | 1     | 0             | 0         | 0         | 0             |  |
| 1        | 0     | 1             | 1         | 1         | 1             | Cannot specify hfname if compact source |
| 1        | 0     | 1             | 1         | 0         | 0             |  |
| 1        | 0     | 1             | 0         | 1         | 1             | Cannot specify hfname if compact source |
| 1        | 0     | 1             | 0         | 0         | 0             |  |
| 1        | 0     | 0             | 1         | 1         | 1             | Cannot have not usecd & not compactsource |
| 1        | 0     | 0             | 1         | 0         | 1             | Cannot have not usecd & not compactsource |
| 1        | 0     | 0             | 0         | 1         | 1             | Cannot have not usecd & not compactsource |
| 1        | 0     | 0             | 0         | 0         | 1             | Cannot have not usecd & not compactsource |
| 0        | 1     | 1             | 1         | 1         | 1             | Cannot specify hfname if compact source |
| 0        | 1     | 1             | 1         | 0         | 0             |  |
| 0        | 1     | 1             | 0         | 1         | 1             | Cannot specify hfname if compact source |
| 0        | 1     | 1             | 0         | 0         | 1             | Must specify efname if compact source and TD-field not specified |
| 0        | 1     | 0             | 1         | 1         | 0             |  |
| 0        | 1     | 0             | 1         | 0         | 1             | If not TD-field and usecd, must specify both efname and hfname |
| 0        | 1     | 0             | 0         | 1         | 1             | If not TD-field and usecd, must specify both efname and hfname |
| 0        | 1     | 0             | 0         | 0         | 1             | If not TD-field and usecd, must specify both efname and hfname |
| 0        | 0     | 1             | 1         | 1         | 1             | Cannot specify hfname if compact source |
| 0        | 0     | 1             | 1         | 0         | 0             |  |
| 0        | 0     | 1             | 0         | 1         | 1             | Cannot specify hfname if compact source |
| 0        | 0     | 1             | 0         | 0         | 1             | Must specify efname if compact source and TD-field not specified |
| 0        | 0     | 0             | 1         | 1         | 1             | Cannot have not usecd & not compactsource |
| 0        | 0     | 0             | 1         | 0         | 1             | Cannot have not usecd & not compactsource |
| 0        | 0     | 0             | 0         | 1         | 1             | Cannot have not usecd & not compactsource |
| 0        | 0     | 0             | 0         | 0         | 1             | Cannot have not usecd & not compactsource |
