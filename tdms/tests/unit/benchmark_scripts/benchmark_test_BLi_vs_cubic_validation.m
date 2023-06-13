%{
BENCHMARK ERROR VALUES FOR BANDLIMITED INTERPOLATION AGAINST CUBIC
INTERPOLATION

Determines the error in interpolated values when using MATLAB's interp
function and cubic interpolation methods, to be used as the benchmark for
test_BLi_vs_cubic_interpolation.cpp.

The function we use is a Cauchy distribution; f(x) = 1/(10x^2+1), over the
range x\in[-2,2]

Neither cubic nor BLi should be the superior interpolation method for all
datapoint-separation  distances.

Results are printed out for each cell size to stdout.
%}
clear;
close all;

% interp inputs
r = 2;
N = 4;
% the spacing between datapoints, mimicing the Yee cell dimensions
cell_Sizes = [0.25, 0.1, 0.05, 0.01];

% for each trial compute BLi and cubic errors
for trial=1:4
    cellSize = cell_Sizes(trial);
    n_datapts = ceil(4/cellSize); % number of Yee cells
    x_lower = -2.; % lower value of the range over which we shall test

    cell_centres = zeros(1,n_datapts-1); % interpolation positions (Yee cell centres)
    field_pos = zeros(1,n_datapts); % data sample positions (field component associated to Yee cells)
    for i=1:n_datapts
        if i~=n_datapts
            cell_centres(i) = x_lower + (i-0.5)*cellSize;
        end
        field_pos(i) = x_lower + (i-1)*cellSize;
    end
    field_samples = BLi_vs_cubic_fn(field_pos); % compute data samples
    exact_vals = BLi_vs_cubic_fn(cell_centres); % compute exact values

    BLi_interp_full = interp(field_samples, r, N); % perform BLi interpolation
    BLi_interp = BLi_interp_full(2:2:end-1); % final point is beyond last Yee cell

    % perform cubic interpolation
    cubic_interp = zeros(1,n_datapts-1);
    cubic_interp(1) = interp_cubic(field_samples(1:4),-1);
    for i=2:n_datapts-2
        cubic_interp(i) = interp_cubic(field_samples(i-1:i+2),0);
    end
    cubic_interp(n_datapts-1) = interp_cubic(field_samples(end-3:end),1);

    % compute errors
    BLi_err = abs( BLi_interp - exact_vals );
    BLi_err_max = max(BLi_err);
    BLi_err_norm = norm( BLi_interp - exact_vals );
    cub_err = abs( cubic_interp - exact_vals );
    cub_err_max = max(cub_err);
    cub_err_norm = norm( cubic_interp - exact_vals );

    BLi_was_better = ones(size(BLi_err));
    BLi_was_better( cub_err < BLi_err ) = 0.;

    fprintf(" ====== \n Cellsize %.3f : \n", cellSize);
    fprintf("MAX Errors: \n");
    fprintf("Cubic: \t %.8e\n", cub_err_max);
    fprintf("BLi: \t %.8e\n", BLi_err_max);
    fprintf("NORM Errors: \n");
    fprintf("Cubic: \t %.8e\n", cub_err_norm);
    fprintf("BLi: \t %.8e\n", BLi_err_norm);
    fprintf("Pointwise, BLi was better at %d points.\n", sum(BLi_was_better));
end

function [out] = interp_cubic(v, pos)
    if pos==0
        % midpoint interpolate
        out = -(1/16) * v(1) + (9/16)*v(2) + (9/16)*v(3) - (1/16)*v(4);
    elseif pos==-1
        % left interp
        out = (5/16)*v(1) + (15/16)*v(2) - (5/16)*v(3) + (1/16)*v(4);
    elseif pos==1
        % right interp
        out = (1/16)*v(1) - (5/16)*v(2) + (15/16)*v(3) + (5/16)*v(4);
    end
end

function [out] = BLi_vs_cubic_fn(x)
    out = 1 ./ (10*x.^2 + 1);
end
