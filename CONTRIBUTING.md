# Contribution Guidelines

### New features/Issues

To suggest a new feature or report an issue please use the 
[issues](https://github.com/UCL/TDMS/issues) board providing details of either 
the feature or the bug.

For bugs, please provide detailed steps to reproduce, your operating system, and what you expected to happen.

To contribute to the code base please
1. [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo)
this repository,
2. commit and push the changes to your fork
3. and submit a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) (PR) against the _main_ branch, being mindful of the guidelines below.

### Pull requests

PRs should follow the following guidelines

- Contain a **context** or **description** of the change
- Include any useful hints for reviewers
- Be focused on a single feature or fix, i.e. a logically consistent set of changes
- Not be excessively large, ideally touching 100s not 1000s of lines
- Add additional tests to ensure the correct functionality
- Should pass all system and unit tests run on GitHub Actions CI
- Follow the code style and have documentation as appropriate

### Tests

Tests are located in [`tdms/tests`](./tdms/tests) and include both unit and system tests. These
must pass locally and on CI for a change to be merged.

### Developer documentation

We have some more detailed developer documentation (C++ class structure, code style, etc) over [here](https://github-pages.ucl.ac.uk/TDMS) on our `gh-pages` site.
