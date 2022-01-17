%function [Cmaterial, Dmaterial, fdtdgrid] = updateupdateterms(Dxl,Dyl, Dzl, dt, delta, material_matrix, composition_matrix, fdtdgrid)
%
%This function sets up the Cmaterial and Dmaterial
%structures. These are structures of the form:
%
%Cmaterial.Ca[x,y,z] - vector length N
%Cmaterial.Cb[x,y,z] - vector length N
%Cmaterial.Cc[x,y,z] - vector length N, only non-zero when we have
%non-dispersive materials
%Dmaterial.Db[x,y,z] - vector length N
%Dmaterial.Da[x,y,z] - vector length N
%
%Where N is the number of materials in the fdtd grid. The material
%field in fdtdgrid will be set to a non zero positive integer if
%the cell is not simply a free space or pml type material. When
%this is the case, the integer is an index into the Cmaterial and
%Dmaterial structures.
%
%see help read_material_data for the specification of material_matrix and composition_matrix
%
%Altered 11/2/2004 to use exponential time stepping for non-zero
%conductivity, this is still experimental
function [Cmaterial, Dmaterial, fdtdgrid] = updateupdateterms(Dxl,Dyl, Dzl, dt, delta, material_matrix, composition_matrix, fdtdgrid)

%define some constants
    [epso , muo , c] = import_constants;

    %store this for later because in my infinite wisdom i over-write eps!
    %eps_0 = eps;
    %matlab now prevents you from using eps both as a function and variable
    eps_0 = 2.2204e-016;

    %The non-uniform cells are specified in a coordinate system with the following origin in global coords
    origin = [Dxl Dyl Dzl];

    %m is the number of cells to alter
    [m,n] = size(composition_matrix);
    [I,J,K] = size(fdtdgrid.Exy);
    [o,p] = size(material_matrix);

    %set up the material structures
    Cmaterial.Cax = zeros(1,o);
    Cmaterial.Cay = zeros(1,o);
    Cmaterial.Caz = zeros(1,o);
    Cmaterial.Cbx = zeros(1,o);
    Cmaterial.Cby = zeros(1,o);
    Cmaterial.Cbz = zeros(1,o);
    Cmaterial.Ccx = zeros(1,o);
    Cmaterial.Ccy = zeros(1,o);
    Cmaterial.Ccz = zeros(1,o);

    Dmaterial.Dax = zeros(1,o);
    Dmaterial.Day = zeros(1,o);
    Dmaterial.Daz = zeros(1,o);
    Dmaterial.Dbx = zeros(1,o);
    Dmaterial.Dby = zeros(1,o);
    Dmaterial.Dbz = zeros(1,o);
    
    if(isempty(composition_matrix) & ~isempty(material_matrix) )
	fprintf(1,'\n\nWarning: An non-empty material matrix has been specified \nalong with an empty composition matrix - the material\nmatrix will not be integrated into the FDTD input file.\n\n');
    end

    if(~isempty(composition_matrix))
	%first check any cells which are out of range
	cellindex = ones(size(composition_matrix,1),1)*origin + composition_matrix(:,1:3);
	criterion = cellindex(:,1) > I | cellindex(:,2) > J | cellindex(:,3) > K;
	inds_out = find(criterion);
	inds_in = find(~criterion);
	for lvar = 1:length(inds_out)
	    fprintf(1,'\nCell specified outside grid - ignoring [%d %d %d]\n',composition_matrix(lvar,1),composition_matrix(lvar,2),composition_matrix(lvar,3));
	end

	%now find unique material identifiers
	[comp_vals,comp_inds] = unique(composition_matrix(:,4));
	identifier_mapping = zeros(max(comp_vals),1);
	identifier_mapping(comp_vals) = 1:length(comp_vals);

	for lvar=1:length(comp_inds)
	    materialindex = comp_vals(lvar);
	    materialprops = material_matrix(max(find(materialindex == material_matrix(:,1))),2:11);
	    materialidentifier = lvar;
	    %%%%
	    eps =  epso*materialprops(1);
	    mu =   muo *materialprops(2);
	    vc =   materialprops(3);
	    wp =   materialprops(4);
	    sigma_x.E = materialprops(5);
	    sigma_y.E = materialprops(6);
	    sigma_z.E = materialprops(7);
	    sigma_x.H = materialprops(8);
	    sigma_y.H = materialprops(9);
	    sigma_z.H = materialprops(10);
	    
	    if max(abs(diff([sigma_x.E sigma_y.E sigma_z.E]))) > 1e-15 |  max(abs(diff([sigma_x.H sigma_y.H sigma_z.H]))) > 1e-15
		
		error('Anisotropic sigmas in interior - program not currently designed to handle this');
		
	    end
	    
	    %Added 11/2/2004 - experimental
	    if(sigma_x.E > eps_0 & abs(vc)<eps_0 & abs(wp)<eps_0)
		fprintf(1,'Non-zero conductivity, using exponential time stepping...');
		
		Cmaterial.Cax(materialidentifier) = exp(-1*sigma_x.E*dt/eps);
		Cmaterial.Cbx(materialidentifier) = (1 - exp(-1*sigma_x.E*dt/eps))./(sigma_x.E*delta.x);
		
		Cmaterial.Cay(materialidentifier) = exp(-1*sigma_y.E*dt/eps);
		Cmaterial.Cby(materialidentifier) = (1 - exp(-1*sigma_y.E*dt/eps))./(sigma_y.E*delta.y);
		
		Cmaterial.Caz(materialidentifier) = exp(-1*sigma_z.E*dt/eps);
		Cmaterial.Cbz(materialidentifier) = (1 - exp(-1*sigma_z.E*dt/eps))./(sigma_z.E*delta.z);
		
	    else
		%end of added section
		
		%handle the dispersive case, will revert to usual case
		%when dispersive parameters are set to 0
		
		gamma = 2*epso*wp*wp*dt*dt/(vc*dt+2);
		
		Cmaterial.Cax(materialidentifier) = (2*eps - sigma_x.E*dt)/(2*eps + .5*gamma + sigma_x.E*dt);
		Cmaterial.Cbx(materialidentifier) = 2*dt/delta.x/(2*eps + 0.5*gamma + sigma_x.E*dt);
		Cmaterial.Ccx(materialidentifier) = 0.5*gamma/(2*eps + 0.5*gamma + sigma_x.E*dt);
		
		Cmaterial.Cay(materialidentifier) = (2*eps - sigma_y.E*dt)/(2*eps + .5*gamma + sigma_y.E*dt);
		Cmaterial.Cby(materialidentifier) = 2*dt/delta.y/(2*eps + 0.5*gamma + sigma_y.E*dt);
		Cmaterial.Ccy(materialidentifier) = 0.5*gamma/(2*eps + 0.5*gamma + sigma_y.E*dt);
		
		Cmaterial.Caz(materialidentifier) = (2*eps - sigma_z.E*dt)/(2*eps + .5*gamma + sigma_z.E*dt);
		Cmaterial.Cbz(materialidentifier) = 2*dt/delta.z/(2*eps + 0.5*gamma + sigma_z.E*dt);
		Cmaterial.Ccz(materialidentifier) = 0.5*gamma/(2*eps + 0.5*gamma + sigma_z.E*dt);
		
% $$$ 		Cmaterial.Cax(materialidentifier) = (1-sigma_x.E*dt/(2*eps))./(1+sigma_x.E*dt/(2*eps));
% $$$ 		Cmaterial.Cbx(materialidentifier) = dt/(eps*delta.x)./ (1+sigma_x.E*dt/(2*eps));
% $$$ 		
% $$$ 		Cmaterial.Cay(materialidentifier) = (1-sigma_y.E*dt/(2*eps))./(1+sigma_y.E*dt/(2*eps));
% $$$ 		Cmaterial.Cby(materialidentifier) = dt/(eps*delta.y)./(1+sigma_y.E*dt/(2*eps));
% $$$ 		
% $$$ 		Cmaterial.Caz(materialidentifier) = (1-sigma_z.E*dt/(2*eps))./(1+sigma_z.E*dt/(2*eps));
% $$$ 		Cmaterial.Cbz(materialidentifier) = dt/(eps*delta.z)./(1+sigma_z.E*dt/(2*eps));
		
	    end
	    
	    
	    
	    Dmaterial.Dax(materialidentifier) = (1-sigma_x.H*dt/(2*mu))./(1+sigma_x.H*dt/(2*mu));
	    Dmaterial.Dbx(materialidentifier) = dt/(mu*delta.x)./(1+sigma_x.H*dt/(2*mu));
	    
	    Dmaterial.Day(materialidentifier) = (1-sigma_y.H*dt/(2*mu))./(1+sigma_y.H*dt/(2*mu));
	    Dmaterial.Dby(materialidentifier) = dt/(mu*delta.y)./(1+sigma_y.H*dt/(2*mu));
	    
	    Dmaterial.Daz(materialidentifier) = (1-sigma_z.H*dt/(2*mu))./(1+sigma_z.H*dt/(2*mu));
	    Dmaterial.Dbz(materialidentifier) = dt/(mu*delta.z)./(1+sigma_z.H*dt/(2*mu));
	    
	    
	end
	inds = sub2ind(size(fdtdgrid.materials),cellindex(:,1),cellindex(:,2),cellindex(:,3));
	%now populate fdtdgrid
	%fdtdgrid.materials(cellindex(:,1),cellindex(:,2),cellindex(:,3)) = identifier_mapping(composition_matrix(:,4));
	%length(find(composition_matrix(:,4)==2))
	fdtdgrid.materials(inds) = identifier_mapping(composition_matrix(:,4));
	%%%%
    end
