# Contribution Guidelines

### New features/Issues

To suggest a new feature or report an issue please use the 
[issue](https://github.com/UCL/TDMS/issues) board providing details of either 
the feature or the bug with steps to reproduce it. To contribute to the code
base please
1. [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo)
this repository,
2. commit and push the changes to your fork
3. and submit a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) (PR) against the _main_ branch, being mindful of the guidelines below.


### Pull requests

PRs should follow the following guidelines

- Contain a **context** or **description** of the change
- Include any useful hints for reviewers
- Not be excessively large (500 line diff at a maximum)
- Add additional tests to ensure the correct functionality
- Should pass all system and unit tests run on GitHub Actions CI


### Tests

Tests are located in [`tdms/tests`](./tdms/tests) and include both unit and system tests. These
must pass locally and on CI for a change to be merged.
