%{
     BENCHMARK FOR E-FIELD INTERPOLATION METHODS

Tests the interpolation accuracy for the E-field components, when
interpolated in each direction to the centre of a Yee cell - benchmarking
test_field_interpolation::TEST_CASE("E-field interpolation check").

The E-field used will have identical field components, each of which are
given by the expression:

E_t(tt) = sin(2\pi tt) exp(-tt^2),

and will be tested over the range $[-2,4]^3$.

Errors are computed as the Frobenius norm of the pointwise difference
(between the exact and interpolated values) matrix, and as the maximum
error in the interpolated values along each axial direction (slice-norm
metric).

Results are displayed to stdout for each component; Ex, Ey, Ez.
%}
clear; 
close all;

% interp parameters
r = 2; 
N = 4;
% bad practice but saves us passing dimensions to subfunctions
% with clear and close all, this should be safe
global cellDims;

% (1,2,3) = (Dx, Dy, Dz)
cellDims = [0.25, 0.1, 0.05];
% domain range is [-2,4]
x_lower = -2.; extent_x = 4;
y_lower = -2.; extent_y = 4;
z_lower = -2.; extent_z = 4;
% function handle with these _lower values
GCC = @(x,y,z) GetCellCentre(x,y,z,x_lower,y_lower,z_lower);
FP = @(x,y,z,c) GetFieldPos(x,y,z,x_lower,y_lower,z_lower,c);

% number of Yee cells in each dimension
nX = round(extent_x/cellDims(1));
nY = round(extent_y/cellDims(2));
nZ = round(extent_z/cellDims(3));

% components of the field at the Yee cell centres,
% and at the field-sample positions
Ex_exact = zeros(nX, nY, nZ); Ex_sample = zeros(nX, nY, nZ);
Ey_exact = zeros(nX, nY, nZ); Ey_sample = zeros(nX, nY, nZ);
Ez_exact = zeros(nX, nY, nZ); Ez_sample = zeros(nX, nY, nZ);
% populate exact values and sample datapoints
for i=1:nX
    for j=1:nY
        for k=1:nZ
            centre = GCC(i,j,k);
            fp_ex = FP(i,j,k,"Ex");
            fp_ey = FP(i,j,k,"Ey");
            fp_ez = FP(i,j,k,"Ez");
            Ex_exact(i,j,k) = Ex_field(centre(1), centre(2), centre(3));
            Ey_exact(i,j,k) = Ey_field(centre(1), centre(2), centre(3));
            Ez_exact(i,j,k) = Ez_field(centre(1), centre(2), centre(3));
            Ex_sample(i,j,k) = Ex_field(fp_ex(1), fp_ex(2), fp_ex(3));
            Ey_sample(i,j,k) = Ey_field(fp_ey(1), fp_ey(2), fp_ey(3));
            Ez_sample(i,j,k) = Ez_field(fp_ez(1), fp_ez(2), fp_ez(3));
        end
    end
end
% NOTE: we don't use the field values at cells with index 0
% so we'll only use E{x,y,z}_exact(2:end,2:end,2:end) from here on out

% interpolate Ex in the x direction to find the values at the cell centres
Ex_interp = zeros(nX-1, nY, nZ);
slice_norms = zeros(nY, nZ);
for k=1:nZ
    for j=1:nY
        temp_x = interp(reshape(Ex_sample(:,j,k),1,nX),r,N);
        % final point is beyond final Yee cell centre
        Ex_interp(:,j,k) = temp_x(2:2:end-1);

        % manual check along this axis of the square-norm error
        slice_norms(j,k) = norm(Ex_interp(:,j,k) - Ex_exact(2:end,j,k));
    end
end
Ex_errs = Ex_interp - Ex_exact(2:end,:,:);
Ex_max_slice_norm = max( slice_norms, [], "all" );
Ex_norm_err = norm( Ex_errs, "fro" );

% interpolate Ey in the y direction to find the values at the cell centres
Ey_interp = zeros(nX, nY-1, nZ);
slice_norms = zeros(nX, nZ);
for k=1:nZ
    for i=1:nX
        temp_y = interp(reshape(Ey_sample(i,:,k),1,nY),r,N);
        % final point is beyond final Yee cell centre
        Ey_interp(i,:,k) = temp_y(2:2:end-1);

        % manual check along this axis of the square-norm error
        slice_norms(i,k) = norm(Ey_interp(i,:,k) - Ey_exact(i,2:end,k));
    end
end
Ey_errs = Ey_interp - Ey_exact(:,2:end,:);
Ey_max_slice_norm = max( slice_norms, [], "all" );
Ey_norm_err = norm( Ey_errs, "fro" );

% interpolate Ez in the z direction to find the values at the cell centres
Ez_interp = zeros(nX, nY, nZ-1);
slice_norms = zeros(nX, nY);
for i=1:nX
    for j=1:nY
        temp_z = interp(reshape(Ez_sample(i,j,:),1,nZ),r,N);
        % final point is beyond final Yee cell centre
        Ez_interp(i,j,:) = temp_z(2:2:end-1);

        % manual check along this axis of the square-norm error
        slice_norms(i,j) = norm(reshape(Ez_interp(i,j,:) - Ez_exact(i,j,2:end), 1, nZ-1));
    end
end
Ez_errs = Ez_interp - Ez_exact(:,:,2:end);
Ez_max_slice_norm = max( slice_norms, [], "all" );
Ez_norm_err = norm( Ez_errs, "fro" );

fprintf("Ex Frobenius norm error: \t %.16e \n", Ex_norm_err);
fprintf("Ex max-slice norm error: \t %.16e \n", Ex_max_slice_norm);
fprintf("Ey Frobenius norm error: \t %.16e \n", Ey_norm_err);
fprintf("Ey max-slice norm error: \t %.16e \n", Ey_max_slice_norm);
fprintf("Ez Frobenius norm error: \t %.16e \n", Ez_norm_err);
fprintf("Ez max-slice norm error: \t %.16e \n", Ez_max_slice_norm);

%% E-field component functions
% Field components E_t(tt) = sin(2\pi tt) exp(-tt^2)
function [value] = Ex_field(~,y,~)
    value = sin(2.*pi*y) * exp(-y.^2);
end
function [value] = Ey_field(~,~,z)
    value = sin(2.*pi*z) * exp(-z.^2);
end
function [value] = Ez_field(x,~,~)
    value = sin(2.*pi*x) * exp(-x.^2);
end

%% Computes the coordinates of the centre of a Yee cell
function [centre] = GetCellCentre(i,j,k,x_lower,y_lower,z_lower)
    global cellDims;
    centre = [x_lower + (i-0.5)*cellDims(1), y_lower + (j-0.5)*cellDims(2), z_lower + (k-0.5)*cellDims(3)];
end

%% Computes the position of the field component associated to a Yee cell
function [fieldPos] = GetFieldPos(i,j,k,x_lower,y_lower,z_lower,component)
    global cellDims;
    fieldPos = GetCellCentre(i,j,k,x_lower,y_lower,z_lower);
    if component=="Ex"
        fieldPos(1) = fieldPos(1) + 0.5*cellDims(1);
    elseif component=="Ey"
        fieldPos(2) = fieldPos(2) + 0.5*cellDims(2);
    elseif component=="Ez"
        fieldPos(3) = fieldPos(3) + 0.5*cellDims(3);
    elseif component=="Hx"
        fieldPos(2) = fieldPos(2) + 0.5*cellDims(2);
        fieldPos(3) = fieldPos(3) + 0.5*cellDims(3);
    elseif component=="Hy"
        fieldPos(1) = fieldPos(1) + 0.5*cellDims(1);
        fieldPos(3) = fieldPos(3) + 0.5*cellDims(3);
    elseif component=="Hz"
        fieldPos(1) = fieldPos(1) + 0.5*cellDims(1);
        fieldPos(2) = fieldPos(2) + 0.5*cellDims(2);
    else
        error("Not a recognised component")
    end
end
