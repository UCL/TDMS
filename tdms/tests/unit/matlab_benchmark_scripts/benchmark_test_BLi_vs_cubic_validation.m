%{
BENCHMARK ERROR VALUES FOR BAND-LIMITED INTERPOLATION AGAINST CUBIC
INTERPOLATION

Determines the error in interpolated values when using MATLAB's interp
function and cubic interpolation methods, to be used as the benchmark for
test_BLi_vs_cubic_interpolation.cpp.

The function we use is a Cauchy distribution multiplied by a sin wave;
f(x) = sin(2\pi x) / (10x^2+1)
over the range x\in[-4,4]

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
cell_Sizes = [7.5e-1, 6.25e-1, 5e-1, 3.75e-1, 2.5e-1, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4];
x_lower = -4.; % lower value of the range over which we shall test
x_upper = 4.; % upper value of the range over which we shall test
extent = x_upper - x_lower;
% make interpolation plots
make_plots = 0;

% for each trial compute BLi and cubic errors
fprintf("Cellsize | (MAX error) BLi :     Cubic       |");
% fprintf(" (Norm error) BLi :     Cubic      |");
fprintf("   (RMSD) BLi    :      Cubic     |");
fprintf(" %%pts BLi better |\n");
fprintf("|:-:|:-:|:-:|:-:|\n");
for trial=1:numel(cell_Sizes)
    cellSize = cell_Sizes(trial);
    n_datapts = ceil(extent/cellSize); % number of Yee cells

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
    BLi_err = BLi_interp - exact_vals;
    BLi_err_max = max(BLi_err);
    BLi_err_norm = norm( BLi_err );
    BLi_rmsd = rmsd(BLi_interp, exact_vals);
    cub_err = cubic_interp - exact_vals;
    cub_err_max = max(cub_err);
    cub_err_norm = norm(cub_err);
    cub_rmsd = rmsd(cubic_interp, exact_vals);

    BLi_was_better = ones(size(BLi_err));
    BLi_was_better( abs(cub_err) < abs(BLi_err) ) = 0.;

    fprintf("%.2e  |", cellSize);
    fprintf("  %.8e : %.8e  |", BLi_err_max, cub_err_max);
%    fprintf("  %.8e  : %.8e |", BLi_err_norm, cub_err_norm);
    fprintf("  %.8e : %.8e |", BLi_rmsd, cub_rmsd);
    fprintf(" %.2f%% |\n", 100. * sum(BLi_was_better) / (n_datapts-1));

    if make_plots
        figure; hold on;
        title(strcat("Interp values: ",num2str(cellSize)));
        plot(cell_centres, BLi_interp);
        plot(cell_centres, cubic_interp);
        plot(cell_centres, exact_vals);
        legend("BLi", "Cubic", "Exact");

        figure; hold on;
        title(strcat("Pt-wise errors: ",num2str(cellSize)));
        plot(cell_centres, BLi_err);
        plot(cell_centres, cub_err);
        legend("BLi", "Cubic");
    end

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
    out = sin(2*pi*x) ./ (10*(x.^2) + 1);
end

% Relative mean square error
function [out] = rmsd(a, b)
    mean_sq_a = mean(abs(a).^2);
    mean_sq_b = mean(abs(b).^2);
    if (mean_sq_a < 1e-16 && mean_sq_b < 1e-16)
        out = 0.;
    else
        out = mean(abs(a-b).^2) / max(mean_sq_a, mean_sq_b);
    end
end
