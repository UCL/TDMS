import matlab.unittest.TestRunner;
import matlab.unittest.Verbosity;
import matlab.unittest.constraints.ContainsSubstring;
import matlab.unittest.selectors.HasBaseFolder;

addpath(genpath('.'));

suite = testsuite(pwd, 'IncludeSubfolders', true);
suite = suite.selectIf(HasBaseFolder(ContainsSubstring('tests')));

runner = TestRunner.withTextOutput('OutputDetail', Verbosity.Detailed);
results = runner.run(suite);
disp(results);

nfailed = nnz([results.Failed]);
assert(nfailed == 0, [num2str(nfailed) ' test(s) failed.']);
