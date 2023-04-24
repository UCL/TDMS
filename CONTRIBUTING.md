# Contribution guidelines

## Bugs

🐛 If you've spotted a bug or want to request a feature please [send an issue through github](https://github.com/UCL/TDMS/issues/new/choose).
Please choose the "bug" or "feature" template (as appropriate) and fill it out with much information as you can.
We really appreciate bug reports and feature requests.
Thanks!

## Code contribution

🚀 If you want to contribute code, or you have reported a bug and think you have a fix, then please:

1. [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo)
this repository,
2. make your changes and check they compile,
3. commit and push the changes to your fork
   + (optional but appreciated) [setup and run pre-commit](https://github-pages.ucl.ac.uk/TDMS/md_doc_developers.html#pre-commit) which automates our code tidyup,
4. and submit a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) (PR) against the `main` branch.

### Pull requests

When you're submitting a PR please give as much context as possible and link to any issues of relevance.
You can use our PR template as a starting point, but here are some guidelines we try to follow:

- Focus on a single feature or bug fix. One set of logically consistent changes per PR.
- Ideally a PR should touch 100s, not 1000s, of lines.
- [Add test(s)](https://github-pages.ucl.ac.uk/TDMS/md_doc_developers.html#testing) to ensure the correct functionality or which reproduce the bug and demonstrate the fix.
- Document your code changes, [we use doxygen for C++](https://github-pages.ucl.ac.uk/TDMS/md_doc_developers.html#code-style-and-doxygen).

You can find more information in our [developer and API documentation](https://github-pages.ucl.ac.uk/TDMS/md_doc_developers.html).
