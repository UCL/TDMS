%function [fdtdgrid] = initialisesplitgrid(I,J,K,Dx,Dy,Dz)
%
%initialises a new FDTD grid for the split formulation with a PML.
%EXCLUDING the PML there are I Yee cells
%in the x-direction, J Yee cells in the y-direction
%and K Yee cells in the z direction. When the PML is taken
%into account there are:
%
%I+2*Dx cells in x direction
%J+2*Dy cells in y direction
%K+2*Dz cells in z direction
%
%Inputs -
%
%I - the number of Yee cells (excluding PML) in the x direction
%J - the number of Yee cells (excluding PML) in the y direction
%K - the number of Yee cells (excluding PML) in the z direction
% Dx - Number of cells in x direction in a single PML layer
% Dy - Number of cells in y direction in a single PML layer
% Dz - Number of cells in z direction in a single PML layer
%
%Returns -
%
% fdtdgrid - a struct with 13 elements:
%               fdtdgrid.Exy
%               fdtdgrid.Exz
%               fdtdgrid.Eyx
%               fdtdgrid.Eyz
%               fdtdgrid.Ezx
%               fdtdgrid.Ezy
%               fdtdgrid.Hxy
%               fdtdgrid.Hxz
%               fdtdgrid.Hyx
%               fdtdgrid.Hyz
%               fdtdgrid.Hzx
%               fdtdgrid.Hzy
%               fdtdgrid.material
%
%These are arrays of diemension (I+2*Dx+1)x(J+2*Dy+1)x(K+2*Dz+1)
%and are intialised to 0. The first 12 are actual field quantities
%whilst the final is a descriptor for a material type. A material
%type of 0 means that the material is either free space or a pml
%and the cell index may be used to index into the appropriate
%vector for the update parameter. If it is greater than zero, the
%update paramter is found by indexing into the material matrix.
%
function [fdtdgrid] = initialisesplitgrid(I,J,K,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu)

I = I + Dxl + Dxu;
J = J + Dyl + Dyu;
K = K + Dzl + Dzu;

fdtdgrid.Exy = zeros(I+1,J+1,K+1);
%fdtdgrid.Exz = zeros(I+1,J+1,K+1);
%fdtdgrid.Eyx = zeros(I+1,J+1,K+1);
%fdtdgrid.Eyz = zeros(I+1,J+1,K+1);
%fdtdgrid.Ezx = zeros(I+1,J+1,K+1);
%fdtdgrid.Ezy = zeros(I+1,J+1,K+1);
%fdtdgrid.Hxy = zeros(I+1,J+1,K+1);
%fdtdgrid.Hxz = zeros(I+1,J+1,K+1);
%fdtdgrid.Hyx = zeros(I+1,J+1,K+1);
%fdtdgrid.Hyz = zeros(I+1,J+1,K+1);
%fdtdgrid.Hzx = zeros(I+1,J+1,K+1);
%fdtdgrid.Hzy = zeros(I+1,J+1,K+1);
fdtdgrid.materials = uint8(zeros(I+1,J+1,K+1));
