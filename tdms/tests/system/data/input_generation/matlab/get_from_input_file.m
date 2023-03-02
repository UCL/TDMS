function [varargout] = get_from_input_file(input_file, optionals, varargin)
    %% Pull the requested variables out of the input_file and return them in the order they were requested.
    % input_file    : Config file to fetch values from
    % optionals     : 1-by-1 struct. The fields of this struct should match the variable names that are considered optional retrievals from the input file, and the values of these fields should be the default value to assign if the variable is not found in the input file.
    %
    % varargin      : Sequence of strings/chars, the names of the variables to retrieve from the input file. They will be returned in the same order as provided to the function as inputs.

%% Allocate storage and useful intermediate values
n_variables = numel(varargin);
varargout = cell(1, n_variables);

optional_vars = fieldnames(optionals);
n_optionals = numel(optional_vars);

%% Load input file
% Check that the input file can be found on the path
if isfile(input_file)
	% Run input file as a script to import variables into the workspace
	run(input_file);
else
	% Throw error - config file cannot be found
    error(sprintf('File %s could not be opened for reading',input_file));
end

%% Determine if the variables loaded in can be retrieved
for i=1:n_variables
    % if the variable was found, save it to the output list
    if exist(varargin{i}, 'var')
        varargout{i} = eval(varargin{i});
    else
        % variable not found, but was it optional?
        for j=1:n_optionals
            if strcmp(varargin{i}, optional_vars{j})
                % this was optional, assign the default value
                varargout{i} = optionals.(optional_vars{j});
                % stop searching to see if this variable was optional
                break
            end
        end
        % we we got to here, the variable was not optional but wasn't found in the input file - throw error
        error(sprintf('Required variable %s not present in %s', required_variables{i}, input_file));
    end
end
end
