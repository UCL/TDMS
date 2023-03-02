function [fdtdgrid] = initialisesplitgrid(I,J,K,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu)
    %% Initialise a new FDTD grid for the split formulation with a PML.
    %% EXCLUDING the PML there are (I,J,K) Yee cells in the (x,y,z)-directions. When the PML is taken into account there are
    %% (I,J,K)+2*(Dx,Dy,Dz)
    %% cells in the (x,y,z) direction.
    % I, J, K       : Number of non- PML Yee cells x, y, z direction
    % Dx, Dy, Dz    : Number of cells in x, y, z direction in a single PML layer
    %
    % fdtdgrid : Struct with 13 elements:
    % Exy, Exz, Eyx, Eyz, Ezx, Ezy, Hxy, Hxz, Hyx, Hyz, Hzx, Hzy : The split-field components
    % material                                                   : Flags whether the cell is PML or not
    %
    % These are arrays of diemension (I+2*Dx+1,J+2*Dy+1,K+2*Dz+1), and are intialised to 0.
    % A material type of 0 means that the material is either free space or a pml cell, and the cell index may be used to index into the appropriate vector for the update parameter.
    % If the material type is greater than zero, the update paramter is found by indexing into the material matrix.

I_tot = I + Dxl + Dxu;
J_tot = J + Dyl + Dyu;
K_tot = K + Dzl + Dzu;

fdtdgrid.Exy = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Exz = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Eyx = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Eyz = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Ezx = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Ezy = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Hxy = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Hxz = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Hyx = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Hyz = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Hzx = zeros(I_tot+1,J_tot+1,K_tot+1);
%fdtdgrid.Hzy = zeros(I_tot+1,J_tot+1,K_tot+1);
fdtdgrid.materials = uint8(zeros(I_tot+1,J_tot+1,K_tot+1));
end
