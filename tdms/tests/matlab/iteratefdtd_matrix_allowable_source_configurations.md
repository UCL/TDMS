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
