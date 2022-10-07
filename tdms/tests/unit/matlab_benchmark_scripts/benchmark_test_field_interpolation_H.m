%{
        BENCHMARK FOR H-FIELD INTERPOLATION METHODS

Tests the interpolation accuracy for the H-field components, when
interpolated in each direction to the centre of a Yee cell - benchmarking
test_field_interpolation::TEST_CASE("H-field interpolation check").

The H-field used will have identical field components, each of which are
given by the expression:

H{x,y,z}(xx,yy,zz) = cos(0.5\pi zz) * exp(-yy^2) * ( 1./ (5xx^2+1) ),

and will be tested over the range $[-2,4]^3$.

Errors are computed as the Frobenius norm of the pointwise difference
(between the exact and interpolated values) matrix. Since we need to slice
in two dimensions for the H-field components, we do not test the slice-norm
metric.

Results are displayed to stdout for each component; Hx, Hy, Hz.
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
x_lower = -2.; extent_x = 4;
y_lower = -2.; extent_y = 4;
z_lower = -2.; extent_z = 4;
% function handle with these _lower values
GCC = @(x,y,z) GetCellCentre(x,y,z,x_lower,y_lower,z_lower);
FP = @(x,y,z,c) GetFieldPos(x,y,z,x_lower,y_lower,z_lower,c);

% number of Yee cells in each dimension
nX = ceil(extent_x/cellDims(1));
nY = ceil(extent_y/cellDims(2));
nZ = ceil(extent_z/cellDims(3));

% components of the field at the Yee cell centres,
% and at the field-sample positions
Hx_exact = zeros(nX, nY, nZ); Hx_sample = zeros(nX, nY, nZ);
Hy_exact = zeros(nX, nY, nZ); Hy_sample = zeros(nX, nY, nZ);
Hz_exact = zeros(nX, nY, nZ); Hz_sample = zeros(nX, nY, nZ);
% populate exact values and sample datapoints
for i=1:nX
    for j=1:nY
        for k=1:nZ
            centre = GCC(i,j,k);
            fp_ex = FP(i,j,k,"Hx");
            fp_ey = FP(i,j,k,"Hy");
            fp_ez = FP(i,j,k,"Hz");
            Hx_exact(i,j,k) = H_field(centre(1), centre(2), centre(3));
            Hy_exact(i,j,k) = H_field(centre(1), centre(2), centre(3));
            Hz_exact(i,j,k) = H_field(centre(1), centre(2), centre(3));
            Hx_sample(i,j,k) = H_field(fp_ex(1), fp_ex(2), fp_ex(3));
            Hy_sample(i,j,k) = H_field(fp_ey(1), fp_ey(2), fp_ey(3));
            Hz_sample(i,j,k) = H_field(fp_ez(1), fp_ez(2), fp_ez(3));
        end
    end
end
% NOTE: we don't use the field values at cells with index 0
% so we'll only use H{x,y,z}_exact(2:end,2:end,2:end) from here on out

% interpolate Hx... do so in y then z, then z then y, and take the worse
% result?
Hx_ypullback = zeros(nX, nY-1, nZ);
Hx_zpullback = zeros(nX, nY, nZ-1);
Hx_interp_yz = zeros(nX, nY-1, nZ-1);
Hx_interp_zy = zeros(nX, nY-1, nZ-1);
for i=1:nX
    for j=1:nY
        temp_z = interp(reshape( Hx_sample(i,j,:),1,nZ),r,N);
        Hx_zpullback(i,j,:) = temp_z(2:2:end-1);
    end
    for k=1:nZ
        temp_y = interp(reshape( Hx_sample(i,:,k),1,nY),r,N);
        Hx_ypullback(i,:,k) = temp_y(2:2:end-1);
    end
    for j=1:nY-1
        temp_z = interp(reshape( Hx_ypullback(i,j,:),1,nZ),r,N);
        Hx_interp_yz(i,j,:) = temp_z(2:2:end-1);
    end
    for k=1:nZ-1
        temp_y = interp(reshape( Hx_zpullback(i,:,k),1,nY),r,N);
        Hx_interp_zy(i,:,k) = temp_y(2:2:end-1);
    end
end
% take worst-case scenario, in the C++ scheme chooses to use one dimension
% over another arbitrarily
Hx_errors = max(abs(Hx_interp_zy - Hx_exact(:,2:end,2:end)), abs(Hx_interp_yz - Hx_exact(:,2:end,2:end)));
Hx_fro_err = norm( Hx_errors, "fro" );
fprintf("Hx Frobenius norm error: \t %.16e \n", Hx_fro_err);

% interpolate Hy... do so in x then z, then z then x, and take the worse
% result?
Hy_xpullback = zeros(nX-1, nY, nZ);
Hy_zpullback = zeros(nX, nY, nZ-1);
Hy_interp_xz = zeros(nX-1, nY, nZ-1);
Hy_interp_zx = zeros(nX-1, nY, nZ-1);
for j=1:nY
    for i=1:nX
        temp_z = interp(reshape( Hy_sample(i,j,:),1,nZ),r,N);
        Hy_zpullback(i,j,:) = temp_z(2:2:end-1);
    end
    for k=1:nZ
        temp_x = interp(reshape( Hy_sample(:,j,k),1,nX),r,N);
        Hy_xpullback(:,j,k) = temp_x(2:2:end-1);
    end
    for i=1:nX-1
        temp_z = interp(reshape( Hy_xpullback(i,j,:),1,nZ),r,N);
        Hy_interp_xz(i,j,:) = temp_z(2:2:end-1);
    end
    for k=1:nZ-1
        temp_x = interp(reshape( Hy_zpullback(:,j,k),1,nX),r,N);
        Hy_interp_zx(:,j,k) = temp_x(2:2:end-1);
    end
end
% take worst-case scenario, in the C++ scheme chooses to use one dimension
% over another arbitrarily
Hy_errors = max(abs(Hy_interp_zx - Hy_exact(2:end,:,2:end)), abs(Hy_interp_xz - Hy_exact(2:end,:,2:end)));
Hy_fro_err = norm( Hy_errors, "fro" );
fprintf("Hy Frobenius norm error: \t %.16e \n", Hy_fro_err);

% interpolate Hz... do so in x then y, then y then x, and take the worse
% result?
Hz_xpullback = zeros(nX-1, nY, nZ);
Hz_ypullback = zeros(nX, nY-1, nZ);
Hz_interp_xy = zeros(nX-1, nY-1, nZ);
Hz_interp_yx = zeros(nX-1, nY-1, nZ);
for k=1:nZ
    for j=1:nY
        temp_x = interp(reshape( Hz_sample(:,j,k),1,nX),r,N);
        Hz_xpullback(:,j,k) = temp_x(2:2:end-1);
    end
    for i=1:nX-1
        temp_y = interp(reshape( Hz_xpullback(i,:,k),1,nY),r,N);
        Hz_interp_xy(i,:,k) = temp_y(2:2:end-1);
    end
    for i=1:nX
        temp_y = interp(reshape( Hz_sample(i,:,k),1,nY),r,N);
        Hz_ypullback(i,:,k) = temp_y(2:2:end-1);
    end
    for j=1:nY-1
        temp_x = interp(reshape( Hz_ypullback(:,j,k),1,nX),r,N);
        Hz_interp_yx(:,j,k) = temp_x(2:2:end-1);
    end
end
Hz_errors = max(abs(Hz_interp_yx - Hz_exact(2:end,2:end,:)), abs(Hz_interp_xy - Hz_exact(2:end,2:end,:)));
Hz_fro_err = norm( Hz_errors, "fro" );
fprintf("Hz Frobenius norm error: \t %.16e \n", Hz_fro_err);

%% H-field component function
% H{x,y,z}(xx,yy,zz) = cos(0.5\pi zz) * exp(-yy^2) * ( 1./ (5xx^2+1) )
function [value] = H_field(x,y,z)
    value = (1 ./ (5 * x.^2 + 1)) * exp(-y.^2) * cos(0.5*pi*z);
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
