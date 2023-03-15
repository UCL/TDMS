%function [fdtdgrid, Ex_out, Ey_out, Ez_out, Hx_out, Hy_out, Hz_out, Ex_bs, Ey_bs, Hx_bs, Hy_bs, x_out, y_out, z_out, Ex_i, Ey_i, Ez_i, Hx_i, Hy_i, Hz_i,x_i,y_i,z_i,vertices,camplitudes,facets,maxresfield] = iteratefdtd_matrix(input_file,operation,outfile,material_file,ill_file)
%
%input_file - file with input configuration information
%operation  - either 'run', 'filesetup', 'gridsetup' or 'illsetup':
%
%       'run': FDTD simulation is directly executed (currently not
%       supported)
%
%       'filesetup': A mat file is written to file outfile, ready
%       for execution.
%
%       'gridsetup': Only the FDTD grid file is setup
%
%       'illsetup': The illumination matrices for pulsed
%       illumination are calculated and saved in matfile given by
%       outfile.
%
%outfile - the filename used to store output, as determined by the
% particular operation chosen.
%
%material_file - the material may be specified as an argument or in
%the input file. If it is specified by both, the one passed as a
%function argument is used.

%ill_file - mat file containing the source field for either pulsed
%illumination or time-domain illumination. If using time-domain
%illumination, this file should contain a struct with two members,
%exi and eyi. These can be either empty or have dimension
%(I+Dxl+Dxu+1) x (J+Dyl+Dyu+1) x Nt
function [fdtdgrid, Ex_out, Ey_out, Ez_out, Hx_out, Hy_out, Hz_out, Ex_bs, Ey_bs, Hx_bs, Hy_bs, x_out, y_out, z_out, Ex_i, Ey_i, Ez_i, Hx_i, Hy_i, Hz_i,x_i,y_i,z_i,vertices,camplitudes,facets,maxresfield] = iteratefdtd_matrix(input_file,operation,outfile,material_file,ill_file)
%Edited 24/7/2003 to_l allow for different cell widths in each orthogonal direction

%% Fetch variables from input file
% Required variables
[delta, I, J, K, n, R0, ...
Dxu, Dxl, Dyu, Dyl, Dzu, Dzl, ...
dt, epsr, mur, f_an, Nt, ...
interface, efname, hfname, wavelengthwidth, ...
illorigin, sourcemode, exphasorsvolume, exphasorssurface, ...
phasorsurface, exdetintegral, outputs_array] = get_from_input_file(input_file, struct(), ...
'delta','I','J','K','n','R0',...
'Dxl','Dxu','Dyl','Dyu','Dzl','Dzu',...
'dt','epsr','mur','f_an','Nt',...
'interface','efname','hfname','wavelengthwidth',...
'illorigin','sourcemode','exphasorsvolume','exphasorssurface',...
'phasorsurface', 'exdetintegral', 'outputs_array');

% Fetch additional requirements if exdetintegral is present and 1
if exdetintegral == 1
	[k_det_obs, NA_det, beta_det, detmodevec, detsensefun] = get_from_input_file(input_file, struct(), 'k_det_obs', 'NA_det', 'beta_det', 'detmodevec', 'detsensefun');
end

% If the user has not specified the matieral file as an input, attempt to obtain it from the input_file
if isempty(material_file)
	material_file = get_from_input_file(input_file, struct(), 'material_file');
	% If there still isn't a non-empty material file, clear the variable name
	if isempty(material_file)
		clear material_file;
	end
end

% Optional variables
% Create struct with default values to assign if not present
optional_args = struct(...
'z_launch', 0, ...
'runmode', 'complete', ...
'dimension', '3', ...
'phasorinc', [1 1 1], ...
'multilayer', [], ...
'kappa_max', 1, ...
'vc_vec', zeros(1, length(epsr)), ...
'wp_vec', zeros(1, length(epsr)), ...
'structure', [], ...
'f_ex_vec', f_an, ...
'intphasorssurface', 1, ...
'exdetintegral', 0, ...
'air_interface', [], ...
'intmatprops', 1, ...
'use_bli', 0, ...
'tdfdir', '', ...
'fieldsample', struct('i',[],'j',[],'k',[],'n',[]), ...
'campssample', struct('vertices',[],'components',[]), ...
'use_pstd', 0 ...
);
[z_launch, runmode, dimension, phasorinc, ...
multilayer, kappa_max, vc_vec, wp_vec, ...
structure, f_ex_vec, intphasorssurface, exdetintegral, ...
air_interface, intmatprops, use_bli, tdfdir, ...
fieldsample, campssample, use_pstd] = get_from_input_file(input_file, optional_args, ...
'z_launch', 'runmode', 'dimension', 'phasorinc', ...
'multilayer', 'kappa_max', 'vc_vec', 'wp_vec', ...
'structure', 'f_ex_vec', 'intphasorssurface', 'exdetintegral', ...
'air_interface', 'intmatprops', 'use_bli', 'tdfdir', ...
'fieldsample', 'campssample', 'use_pstd');

%% Validate variables in the workspace
% Check that the outputs_array is valid
for lvar = 1:size(outputs_array,1)
    if length(outputs_array{lvar}) < 6
        error('insufficient number of entries in outputs_array number %d',lvar);
    end

	% Check range vectors
	node1 = outputs_array{lvar}{3};
	node2 = outputs_array{lvar}{4};
	if ~all(node1(1:3)<=node2(1:3))
		error('outputs_array entry %d range vectors incorrect',lvar);
	end

	% Adjust to global indexing if necessary
	lower_dim_limit = [1 1 1];
	if strcmp(outputs_array{lvar}{2},'interior')
		% Indexing is relative to the first non-pml cell
		if strcmp(dimension, '3')
			upper_dim_limit = [I J K];
		else
			upper_dim_limit = [I J 1];
		end
		% Set values to use global index
		outputs_array{lvar}{3} = outputs_array{lvar}{3} + [Dxl Dyl Dzl];
		outputs_array{lvar}{4} = outputs_array{lvar}{4} + [Dxl Dyl Dzl];
	else
		% Indexing is global, relative to 1st cell in PML grid
		% TSTK - shouldn't we check for dimension 0 here too (K -> 1?)
		upper_dim_limit = [I J K] + [Dxl Dyl Dzl] + [Dxu Dyu Dzu] + 1;
	end
	% Check that cell indices were within the allowed range
	% Compare to node1 and node2 in case we hae already adjusted to global coordinates
	if all(node1 <= upper_dim_limit) && all(node1 >= lower_dim_limit) && all(node2 <= upper_dim_limit) && all(node2 >= lower_boundary_update)
		error('outputs_array entry %d range vectors incorrect',lvar);
	end

	% Check components are valid
	components_to_fetch = outputs_array{lvar}{5};
	if ~isempty(components_to_fetch)
		% If components is non-empty, but also doesn't specify (at least one of) the x, y, z, components, throw error
		if ~(contains(components_to_fetch, 'x') || contains(components_to_fetch, 'y') || contains(components_to_fetch, 'z'))
			error('%s is not a vector component',  components_to_fetch);
		end
	else
		% No components are specified, raise an error
		error('No components were specified in %s (outputs_array{%d}', components_to_fetch, lvar);
	end

	% Check output style is accumulate or dump
	output_style = outputs_array{lvar}{6};
	if ~(strcmp(output_style, 'accumulate') || strcmp(output_style, 'dump'))
		error('Mode of output writing (%s) is invalid for outputs_array{%d}', output_style, lvar);
	end
end

% Check campssample has correct fields
if ~isempty(campssample)
	if numel(fieldnames(campssample)) ~= 2
		error('campssample must have exactly two fields (components and vertices)');
	elseif ~isfield(campssample, 'components')
		error('campssample does not have a components field');
	elseif ~isfield(campssample, 'verticies')
		error('camppsample does not have a vertices field');
	end

	% Cast to expected datatype
    campssample.vertices = int32( campssample.vertices );
    campssample.components = int32( campssample.components );
end

% Check that interface has been fully specified.
% If not, set any un-specified interface conditions 0.
% That is, the interface condition will not be implemented
interface_field = {'I0','I1','J0','J1','K0','K1'};
for lvar = 1:length(interface_field)
    if ~isfield(interface,interface_field{lvar})
		fprintf('interface.%s has not been defined, setting it to [0 0]\n', interface_field{lvar});
		interface.(interface_field{lvar}) = [0 0];
    end
end
% Now check that the entries for I0...K1 are valid
interface_component_valid(interface.I0, interface.I1, I);
interface_component_valid(interface.J0, interface.J1, J);
interface_component_valid(interface.K0, interface.K1, K);

% Check that sourcemode is valid
if ~(strcmp(sourcemode,'pulsed') || strcmp(sourcemode,'steadystate'))
    error('sourcemode (%s) invalid', sourcemode);
end

% Check that operation is valid, and based on this whether runmode is valid
if ~(strcmp(operation,'run') || strcmp(operation,'filesetup') || strcmp(operation,'gridsetup') || strcmp(operation,'illsetup'))
    error('operation (%s) invalid',operation);
elseif strcmp(operation,'filesetup') && ~strcmp(runmode, 'complete')
	error('Only runmode complete is suported for operation filesetup');
else
	% Check that runmode is valid
	if ~(strcmp(runmode,'analyse') || strcmp(runmode,'complete'))
		error('runmode (%s) invalid', runmode);
	end
end

% Check that phasorsurface is valid
if exphasorssurface
	if numel(phasorsurface) ~= 6
		error('phasorsurface must be a vector of 6 elements');
	else
		% Now check that each entry is valid
		% For J, we also need to confirm this isn't a 2D simulation (so there's >1 cell in the J-direction)
		not_2D_simulation = ((J + Dyl + Dyu) > 0);
		% Otherwise, for each component x/y/z and the corresponding number of cells in that dimension I/J/K, we need
		% 1 <= phasorsurface(lower-index) + D{x,y,z}l - 2,
		% phasorsurface(upper-index) - D{x,y,z}u <= {I,J,K}
		if phasorsurface(1) + Dxl -2 < 1
			error('phasorsurface(1) out of range');
		elseif phasorsurface(2) - Dxu > I
			error('phasorsurface(2) out of range');
		elseif (phasorsurface(3) + Dyl -2 < 1) && not_2D_simulation
			error('phasorsurface(3) out of range');
		elseif (phasorsurface(4) - Dyu > J) && not_2D_simulation
			error('phasorsurface(4) out of range');
		elseif (phasorsurface(5) + Dzl -2 < 1) && ~( (K==0) && (phasorsurface(5)==0) )
			error('phasorsurface(5) out of range');
		elseif phasorsurface(6) - Dzu > K
			error('phasorsurface(6) out of range');
		end
		% Now convert to global coordinates (still in MATLAB 1-indexing though)
		phasorsurface = phasorsurface + [Dxl Dxl Dyl Dyl Dzl Dzl];
	end
end
% Now check that the phasor surface has been specified correctly
if ( mod( (phasorsurface(2) - phasorsurface(1)),phasorinc(1) ) || ...
	 mod( (phasorsurface(4) - phasorsurface(3)),phasorinc(2) ) || ...
	 mod( (phasorsurface(6) - phasorsurface(5)),phasorinc(3) ) )
    error('Incorrect specification of phasorinc');
end

% Cast band-limited to a boolean value
use_bli = all(logical(use_bli));
% Cast use_pstd to a boolean value
use_pstd = all(logical(use_pstd));

% Check dimension is valid
if ~(strcmp('3', dimension) || strcmp('TE', dimension) || strcmp('TM', dimension))
    error('Dimension must be either ''3'', ''TM'' or ''TE''');
end

% Compute the extent of the grid - including the PML
I_tot = I + Dxl + Dxu;
J_tot = J + Dyl + Dyu;
K_tot = K + Dzl + Dzu;

% Convert campssample.vertices to global coordinates, MATLAB version
if ~isempty(campssample.vertices)
    campssample.vertices(:,1) = campssample.vertices(:,1) + Dxl;
    campssample.vertices(:,2) = campssample.vertices(:,2) + Dyl;
    campssample.vertices(:,3) = campssample.vertices(:,3) + Dzl;
else
	% If vertices is empty, make sure components is also empty
    campssample.components = [];
end

% Validate structure
if isempty(multilayer) && ~isempty(structure)
	error('structure is not empty, but multlayer is');
elseif ~isempty(structure)
	structure(1,:) = structure(1,:) + Dxl;
	old_structure = structure;
	if old_structure(1,1)>1
		structure = zeros(2,size(structure,2)+1);
		structure(:,2:end) = old_structure;
		structure(1,1) = 1;
		structure(2,1) = old_structure(2,1);
	end

	old_structure = structure;
	if old_structure(1,end)<I_tot
		structure = zeros(2,size(structure,2)+1);
		structure(:,1:(end-1)) = old_structure;
		structure(1,end) = I_tot+1;
		structure(2,end) = old_structure(2,end);
	end

	tmp_pro = cumprod(double(diff(structure(1,:))>0));
	if  ~tmp_pro(end)
		error('structure should have first row monotonically increasing');
	elseif structure(1,1)<1
		error('structure should have first row entries all greater than 1');
	elseif structure(1,end)>(I_tot+1)
		error('structure should have first row entries all less than I+1');
	end

	[ml_mat,pl_mat] = ndgrid(multilayer,structure(2,:));
	if ~isempty(find(ml_mat+pl_mat<2, 1))
		error('trench is breaking PML boundary, reduce magnitude of structure displacement');
	elseif ~isempty(find(ml_mat+pl_mat>K-1, 1))
		error('trench is breaking PML boundary, reduce magnitude of structure displacement');
	end

	i_int = 1:(I_tot+1);
	p_int = int32(round(interp1(structure(1,:),structure(2,:),i_int,'linear')));
	i_int = int32(i_int);

	structure = [i_int;p_int];
end

% Validate multilayer
if ~isempty(multilayer)
    if min(multilayer) < 1
	error('multilayer must have elements greater than 0');
	elseif max(multilayer) > K
		error('multilayer must have elements less than or equal to  K');
    elseif all(sort(multilayer)==multilayer)
		error('multilayer must be monotonically increasing');
    end
end

% Check kappa_max dimension
if numel(kappa_max)~=1 && numel(kappa_max)~=numel(epsr)
	error('should have numel(kappa_ax)=numel(epsr)');
end

% Check vc_vec and wp_vec have identical dimension, and one greater than that of multilayer
if numel(vc_vec) ~= numel(wp_vec)
    error('vc_vec and wp_vec should have the same number of elements');
elseif numel(multilayer)~=(numel(vc_vec)-1)
    error('vc_vec and wp_vec should have one more member than multilayer');
end

% Check that multilayer and epsr have the correct dimension(s)
if isempty(multilayer) && numel(epsr)~=1
	error('epsr should have only a single element in the case of a single layer');
elseif numel(multilayer)~=(numel(epsr)-1)
	error('epsr should have one more member than multilayer');
end

% Check f_ex_vec is either a row or column vector
if numel(f_ex_vec) > 1 && ~(size(f_ex_vec,1)==1 || size(f_ex_vec,2)==1)
	error('f_ex_vec should be a vector (ie, a matrix with one singleton dimension)');
end

%check correctness of detector sensitivity evaluation parameters
if ~(exdetintegral==0 || exdetintegral==1)
    error('exdetintegral should take a boolean value (0 or 1)');
elseif exdetintegral==1
    if NA_det < 0 || NA_det > 1
		error('NA_det should be between 0 and 1');
	elseif sum(abs( round(detmodevec) - detmodevec ))~=0
		error('detmodevec should contain only integers');
    elseif ~all(detmodevec > 0)
		error('detmodevec should contain only positive integers');
	elseif isempty(detmodevec)
		error('detmodevec should contain at least one element');
    end
	% Regardless of shape, cast this to a column vector array
    detmodevec = detmodevec(:);
end

% Check that tdfdir does not end with a (forward)slash
if(~isempty(tdfdir)) && tdfdir(end)=='/'
	tdfdir = tdfdir(1:(end-1));
end

% Check air_interface
if ~isempty(air_interface) && air_interface > -illorigin(3)*delta.z
	error('Error: air_interface is intended to be defined outside of the simulation space (z=%e), in the negative z-direction, however it is at z=%e',-illorigin(3)*delta.z,air_interface);
end

%% End of error checking, begin setup of .mat file output

fprintf('Allocating grid...');
fdtdgrid = initialisesplitgrid(I,J,K,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu);
fprintf('Done\n');

[epso , ~ , c] = import_constants;

% We now implement an incident planewave using a total/scattered field formulation.
% For k<K1 we have scattered and for k>=K1 we have total.
% The field is coupled into system by means of update equations to_l maintain consistency.

refractive_index = sqrt(real(epsr(1)));
omega_an = 2*pi*f_an;
lambda_an = c/(f_an*refractive_index);

fprintf('Initialising source field...\n');

% Offset the illumination origin to be in global coordinates
illorigin = illorigin + [Dxl Dyl Dzl];
% Now need to adjust interface to be in global coordinates, unless we are using a pulsed source.
if strncmp(sourcemode,'pulsed',6)
	% For a pulsed source, the interfaces need to be hard-reset
    interface.I0(1) = 1;
    interface.J0(1) = 1;
    interface.I1(1) = I_tot+1;
    interface.J1(1) = J_tot+1;

    interface.I0(2) = 0;
    interface.I1(2) = 0;
    interface.J0(2) = 0;
    interface.J1(2) = 0;
    interface.K1(2) = 0;

    if (interface.K0(2)==0) && K~=0
	    error('Running in pulsed mode with k0[0]=0, there is no point running');
    end
else
	interface.I0(1) = interface.I0(1) + Dxl;
	interface.I1(1) = interface.I1(1) + Dxl;
	interface.J0(1) = interface.J0(1) + Dyl;
	interface.J1(1) = interface.J1(1) + Dyl;
	interface.K0(1) = interface.K0(1) + Dzl;
	interface.K1(1) = interface.K1(1) + Dzl;
end

% Setup the source terms
if length(ill_file) > 0 %must have already computed the illumination source
    fprintf('Loading illumination source from %s\n', ill_file);

    data = load(ill_file);
    %here we can have a data file with elemenets Isource, Jsource
    %and Ksource *or* exi and eyi
    fieldnames_ill = fieldnames(data);

    if has_ijk_source_matricies(data)
        Isource = data.Isource;
        Jsource = data.Jsource;
        Ksource = data.Ksource;

        tdfield.exi = [];
        tdfield.eyi = [];
        assert_source_has_correct_dimensions(Isource, Jsource, Ksource, interface);
    elseif has_exi_eyi(data)

		non_empty_eh_names = (~isempty(efname)) && (~isempty(hfname));
        tdfield = data;
        assert_exi_eyi_have_correct_dimensions(data.exi, data.eyi, I_tot, J_tot, Nt);

        if (interface.I0(2) | interface.I1(2)) && non_empty_eh_names
            Isource = zeros(8,interface.J1(1) - interface.J0(1) + 1, interface.K1(1) - interface.K0(1) + 1);
        else
            Isource=[];
        end

        if (interface.J0(2) | interface.J1(2)) && non_empty_eh_names
            Jsource = zeros(8,interface.I1(1) - interface.I0(1) + 1, interface.K1(1) - interface.K0(1) + 1);
        else
            Jsource = [];
        end

        if (interface.K0(2) | interface.K1(2)) && non_empty_eh_names
            Ksource = zeros(8,interface.I1(1) - interface.I0(1) + 1, interface.J1(1) - interface.J0(1) + 1);
        else
            Ksource = [];
        end

		if non_empty_eh_names
			%Set up the Isource field. This has to be defined on a 2d array
			%over the range (J0,J1)x(K0,K1). We calculate the field values
			%assuming an origin for the illumination
			i_source = interface.I0(1) - illorigin(1);
			j_source = (interface.J0(1):interface.J1(1)) - illorigin(2);
			k_source = (interface.K0(1):interface.K1(1)) - illorigin(3);

			%Ey, I0
			if interface.I0(2)
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				Isource(1,:,:) = source_field{2};

				%Ez, I0
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ez');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				Isource(2,:,:) = source_field{3};

				%Hy, I0
				[x,y,z] = yeeposition(i_source-1,j_source,k_source,delta,'Hy');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Isource(3,:,:) = source_field{2};

				%Hz, I0
				[x,y,z] = yeeposition(i_source-1,j_source,k_source,delta,'Hz');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Isource(4,:,:) = source_field{3};
			end

			i_source = interface.I1(1) - illorigin(1);
			%Ey, I1
			if interface.I1(2)
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				Isource(5,:,:) = source_field{2};

				%Ez, I1
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ez');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				Isource(6,:,:) = source_field{3};

				%Hy, I1
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hy');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Isource(7,:,:) = source_field{2};

				%Hz, I1
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hz');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Isource(8,:,:) = source_field{3};
			end

			%Set up the Jsource field. This has to be defined on a 2d array
			%over the range (I0,I1)x(K0,K1). We calculate the field values
			%assuming an origin for the illumination
			i_source = (interface.I0(1):interface.I1(1)) - illorigin(1);
			j_source = interface.J0(1) - illorigin(2);
			k_source = (interface.K0(1):interface.K1(1))- illorigin(3);


			%Ey, J0
			if interface.J0(2)
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				Jsource(1,:,:) = source_field{1};

				%Ez, J0
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ez');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				Jsource(2,:,:) = source_field{3};

				%Hy, J0
				[x,y,z] = yeeposition(i_source,j_source-1,k_source,delta,'Hx');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Jsource(3,:,:) = source_field{1};

				%Hz, J0
				[x,y,z] = yeeposition(i_source,j_source-1,k_source,delta,'Hz');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Jsource(4,:,:) = source_field{3};
			end

			j_source = interface.J1(1) - illorigin(2);
			%Ey, J1
			if interface.J1(2)
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				Jsource(5,:,:) = source_field{1};

				%Ez, J1
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ez');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				Jsource(6,:,:) = source_field{3};

				%Hy, J1
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hx');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Jsource(7,:,:) = source_field{1};

				%Hz, J1
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hz');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Jsource(8,:,:) = source_field{3};
			end

			%Set up the Ksource field. This has to be defined on a 2d array
			%over the range (I0,I1)x(J0,J1). We calculate the field values
			%assuming an origin for the illumination
			i_source = (interface.I0(1):interface.I1(1)) - illorigin(1);
			j_source = (interface.J0(1):interface.J1(1)) - illorigin(2);
			k_source = interface.K0(1) - illorigin(3);


			%Ex, K0
			if interface.K0(2)
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
				z = z + z_launch;
				%fprintf(1,'%d %d %d %e\n',i_source,j_source,k_source,z);
				[X,Y,Z] = ndgrid(x,y,z);
				if use_pstd
					eval(sprintf('source_field = %s(X,Y,Z-delta.z/2);',efname ));
				else
					eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				end
				Ksource(1,:,:) = source_field{1};

				%Ey, K0
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				if use_pstd
					eval(sprintf('source_field = %s(X,Y,Z-delta.z/2);',efname ));
				else
					eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				end
				Ksource(2,:,:) = source_field{2};
				%Hx, K0
				[x,y,z] = yeeposition(i_source,j_source,k_source-1,delta,'Hx');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Ksource(3,:,:) = source_field{1};

				%Hy, K0
				[x,y,z] = yeeposition(i_source,j_source,k_source-1,delta,'Hy');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Ksource(4,:,:) = source_field{2};
			end

			k_source = interface.K1(1) - illorigin(3);
			%Ex, K1
			if interface.K1(2)
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				Ksource(5,:,:) = source_field{1};

				%Ey, K1
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
				Ksource(6,:,:) = source_field{2};

				%Hx, K1
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hx');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Ksource(7,:,:) = source_field{1};

				%Hy, K1
				[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hy');
				z = z + z_launch;
				[X,Y,Z] = ndgrid(x,y,z);
				eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
				Ksource(8,:,:) = source_field{2};
			end
		end
	else
        error('TDMSException:InvalidIlluminationFile', ...
             'Illumination file did not have the correct elements. Need either {Isource, Jsoruce, Ksource} or {exi, eyi}');
    end
else
	% an illumination file has not been specified

    fprintf('Creating Isource, Jsource, Ksource...');

	non_empty_eh_names = (~isempty(efname)) && (~isempty(hfname));
    tdfield.exi = [];tdfield.eyi = [];
    if (interface.I0(2) | interface.I1(2)) && non_empty_eh_names
		Isource = zeros(8,interface.J1(1) - interface.J0(1) + 1, interface.K1(1) - interface.K0(1) + 1);
    else
		Isource=[];
    end

    if (interface.J0(2) | interface.J1(2)) && non_empty_eh_names
		Jsource = zeros(8,interface.I1(1) - interface.I0(1) + 1, interface.K1(1) - interface.K0(1) + 1);
    else
		Jsource = [];
    end

    if (interface.K0(2) | interface.K1(2)) && non_empty_eh_names
		Ksource = zeros(8,interface.I1(1) - interface.I0(1) + 1, interface.J1(1) - interface.J0(1) + 1);
    else
		Ksource = [];
    end


    %Set up the Isource field. This has to be defined on a 2d array
    %over the range (J0,J1)x(K0,K1). We calculate the field values
    %assuming an origin for the illumination
    i_source = interface.I0(1) - illorigin(1);
    j_source = (interface.J0(1):interface.J1(1)) - illorigin(2);
    k_source = (interface.K0(1):interface.K1(1)) - illorigin(3);

	if non_empty_eh_names
		%Ey, I0
		if interface.I0(2)
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			Isource(1,:,:) = source_field{2};

			%Ez, I0
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ez');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			Isource(2,:,:) = source_field{3};

			%Hy, I0
			[x,y,z] = yeeposition(i_source-1,j_source,k_source,delta,'Hy');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Isource(3,:,:) = source_field{2};


			%Hz, I0
			[x,y,z] = yeeposition(i_source-1,j_source,k_source,delta,'Hz');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Isource(4,:,:) = source_field{3};
		end

		i_source = interface.I1(1) - illorigin(1);
		%Ey, I1
		if interface.I1(2)
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			Isource(5,:,:) = source_field{2};

			%Ez, I1
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ez');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			Isource(6,:,:) = source_field{3};

			%Hy, I1
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hy');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Isource(7,:,:) = source_field{2};

			%Hz, I1
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hz');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Isource(8,:,:) = source_field{3};
		end

		%Set up the Jsource field. This has to be defined on a 2d array
		%over the range (I0,I1)x(K0,K1). We calculate the field values
		%assuming an origin for the illumination
		i_source = (interface.I0(1):interface.I1(1)) - illorigin(1);
		j_source = interface.J0(1) - illorigin(2);
		k_source = (interface.K0(1):interface.K1(1))- illorigin(3);


		%Ey, J0
		if interface.J0(2)
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			Jsource(1,:,:) = source_field{1};

			%Ez, J0
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ez');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			Jsource(2,:,:) = source_field{3};

			%Hy, J0
			[x,y,z] = yeeposition(i_source,j_source-1,k_source,delta,'Hx');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Jsource(3,:,:) = source_field{1};

			%Hz, J0
			[x,y,z] = yeeposition(i_source,j_source-1,k_source,delta,'Hz');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Jsource(4,:,:) = source_field{3};
			end

		j_source = interface.J1(1) - illorigin(2);
		%Ey, J1
		if interface.J1(2)
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			Jsource(5,:,:) = source_field{1};

			%Ez, J1
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ez');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			Jsource(6,:,:) = source_field{3};

			%Hy, J1
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hx');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Jsource(7,:,:) = source_field{1};

			%Hz, J1
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hz');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Jsource(8,:,:) = source_field{3};
		end

		%Set up the Ksource field. This has to be defined on a 2d array
		%over the range (I0,I1)x(J0,J1). We calculate the field values
		%assuming an origin for the illumination
		i_source = (interface.I0(1):interface.I1(1)) - illorigin(1);
		j_source = (interface.J0(1):interface.J1(1)) - illorigin(2);
		k_source = interface.K0(1) - illorigin(3);


		%Ex, K0
		if interface.K0(2)
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
			z = z + z_launch;
			%fprintf(1,'%d %d %d %e\n',i_source,j_source,k_source,z);
			[X,Y,Z] = ndgrid(x,y,z);
			if use_pstd
				eval(sprintf('source_field = %s(X,Y,Z-delta.z/2);',efname ));
			else
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			end
			Ksource(1,:,:) = source_field{1};

			%Ey, K0
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			if use_pstd
				eval(sprintf('source_field = %s(X,Y,Z-delta.z/2);',efname ));
			else
				eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			end
			Ksource(2,:,:) = source_field{2};
			%Hx, K0
			[x,y,z] = yeeposition(i_source,j_source,k_source-1,delta,'Hx');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Ksource(3,:,:) = source_field{1};

			%Hy, K0
			[x,y,z] = yeeposition(i_source,j_source,k_source-1,delta,'Hy');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Ksource(4,:,:) = source_field{2};
			end

		k_source = interface.K1(1) - illorigin(3);
		%Ex, K1
		if interface.K1(2)
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			Ksource(5,:,:) = source_field{1};

			%Ey, K1
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',efname ));
			Ksource(6,:,:) = source_field{2};

			%Hx, K1
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hx');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Ksource(7,:,:) = source_field{1};

			%Hy, K1
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hy');
			z = z + z_launch;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('source_field = %s(X,Y,Z);',hfname ));
			Ksource(8,:,:) = source_field{2};
		end
	end
end
fprintf('Done initialising source field\n');

if strncmp(operation,'illsetup',8)%save the source terms
    save(outfile,'Isource','Jsource','Ksource');
    Ex_out = [];
    Ey_out = [];
    Ez_out = [];

    Hx_out = [];
    Hy_out = [];
    Hz_out = [];

    Ex_bs = [];
    Ey_bs = [];


    Hx_bs = [];
    Hy_bs = [];

    x_out = [];
    y_out = [];
    z_out = [];

    Ex_i = [];
    Ey_i = [];
    Ez_i = [];
    Hx_i = [];
    Hy_i = [];
    Hz_i = [];
    x_i = [];
    y_i = [];
    z_i = [];

    vertices = [];
    camplitudes = [];
    facets = [];
    maxresfield = [];
else
    if J==0
		dt_upper = 1/(c*sqrt(1/delta.x^2 + 1/delta.z^2))*min(sqrt(real(epsr)));
    else
		dt_upper = 1/(c*sqrt(1/delta.x^2 + 1/delta.y^2 + 1/delta.z^2))*min(sqrt(real(epsr)));
    end

    if dt > dt_upper
		fprintf(1, 'dt being changed from %.5e to %.5e to ensure stability\n',dt,dt_upper*0.95);
		dt = dt_upper*0.95;
    end
    %setup a guassian pulse g(t) = exp(-pi*(( (t-to_l)/(hwhm) )^2))
    hwhm = lambda_an^2/((c/refractive_index)*wavelengthwidth)*2*sqrt(log(2)/pi);
    to_l = hwhm*sqrt(log(1e8)/pi);
    %to_u = to_l + (K1-K0)*delta.z/c;
    to_u = to_l + (interface.K1(1)-interface.K0(1))*delta.z/c;

    %complete_pulse = 1 means that a guassian pulse is launched. complete_pulse = 0
    %means that the leading edge of the incident field is guassian however incident
    %field becomes steady state when the pulse reaches its peak
    complete_pulse = 1;

    %generate the update constants
    fprintf('Initialising update terms [free space and pml]...');

    [sigma, C, D, freespace, conductive_aux, dispersive_aux] = initialiseupdateterms(R0, I, J, K, Dxl, ...
						     Dxu, Dyl, Dyu, Dzl, ...
						     Dzu, n, delta, dt, ...
						     epsr, mur, multilayer, omega_an, kappa_max,vc_vec,wp_vec);
    fprintf('Done\n');

    fprintf('Initialising update terms [reading grid composition from %s]...',material_file);
    iteratefdtd_matrix_path = which('iteratefdtd_matrix');
    file_parser_path = sprintf('%s/file_parser',iteratefdtd_matrix_path(1:(findstr(iteratefdtd_matrix_path,'iteratefdtd_matrix')-1)));
    addpath(file_parser_path);
    [material_matrix,composition_matrix] = read_material_data(material_file);
    rmpath(file_parser_path);

    %now update the update parameters to take account of the material
    %composition specified in material_file

    [Cmaterial, Dmaterial, fdtdgrid] = updateupdateterms(Dxl,Dyl, Dzl, dt, delta, material_matrix, composition_matrix, fdtdgrid);

    %now we setup the alpha, beta and gamma terms for the
    %dispersive materials.

    %first we have to so some sorting and ordering since the
    %material matrix might not be in order and may have gaps
    if ~isempty(material_matrix)
		mat_inds = material_matrix(:,1);

		vc = zeros(1,max(mat_inds));
		wp = zeros(1,max(mat_inds));
		epsr_vec = zeros(1,max(mat_inds));


		vc(mat_inds) = material_matrix(:,4);
		wp(mat_inds) = material_matrix(:,5);
		epsr_vec(mat_inds) = material_matrix(:,2);

		disp_params.alpha = 4./(vc*dt + 2);
		disp_params.beta = (vc*dt-2)./(vc*dt+2);
	%	disp_params.gamma =2.*epsr_vec*epso.*wp.*wp*dt*dt./(vc*dt+2);
			disp_params.gamma =2.*epso.*wp.*wp*dt*dt./(vc*dt+2);
    else
		disp_params.alpha = [];
		disp_params.beta =  [];
		disp_params.gamma = [];
    end

    fprintf('Done\n');

    %time vectors
    %tvec_E = (1:Nt)*dt;     %tvec_E is the time that the incident E-field is required
    %tvec_H = tvec_E-dt/2;   %tvec_H is the time that the incident H-field is required

    %now set up the cartesian axis coordinates of the entire grid
    %including the PML. This is done so that the internal indexing
    %scheme of the FDTD grid also indexes the coordinates of the grid

    x_grid_label = ((1:(I_tot + 1)) - illorigin(1))*delta.x;
    y_grid_label = ((1:(J_tot + 1)) - illorigin(2))*delta.y;
    z_grid_label = ((1:(K_tot + 1)) - illorigin(3))*delta.z + z_launch;
    grid_labels.x_grid_labels = x_grid_label;
    grid_labels.y_grid_labels = y_grid_label;
    grid_labels.z_grid_labels = z_grid_label;

    [m_outputs,n_outputs] = size(outputs_array);
    m_outputs
    Nt
    accumulate_output = cell(m_outputs, Nt);

    %now set the fdtdgrid to an initial state
    fprintf('Initialising the field values...');
    t_init_E = 0;
    t_init_H = 0.5*dt;
    phasetermE = exp(sqrt(-1)*omega_an*t_init_E);
    phasetermH = exp(sqrt(-1)*omega_an*t_init_H);

    fprintf('Done\n');
    %end setting fdtdgrid

    %now setup detector integral data, if required
    if exdetintegral == 1
		%fx_vec, fy_vec
		fx_vec = (0:(I-1))/I/delta.x;
		fy_vec = (0:(J-1))/J/delta.y;

		fx_vec(find(fx_vec>1/2/delta.x)) = fx_vec(find(fx_vec>1/2/delta.x)) - 1/delta.x;
		fy_vec(find(fy_vec>1/2/delta.y)) = fy_vec(find(fy_vec>1/2/delta.y)) - 1/delta.y;

		f_vec.fx_vec = fx_vec;
		f_vec.fy_vec = fy_vec;

		%FX, FY
		[FX,FY]=ndgrid(fx_vec,fy_vec);

		%Pupil - correct for for both cases of with and without an
			%air interface
		Pupil = double( sqrt( (lambda_an*FX).^2+(lambda_an*FY).^2)< NA_det/refractive_index );

		%k_det_obs_global
		k_det_obs_global = k_det_obs + Dzl;

		%Dx_tilde, Dy_tilde
		%Handle 2D case
		if I==0
			i_obs = 1-(illorigin(1)-Dxl);
		else
			i_obs = (1:I)-(illorigin(1)-Dxl);
		end

		if J==0
			j_obs = 1-(illorigin(2)-Dyl);
		else
			j_obs = (1:J)-(illorigin(2)-Dyl);
		end

		for im=1:numel(detmodevec)

			[x,y,z] = yeeposition(i_obs,j_obs,k_det_obs-illorigin(3)+Dzl,delta,'Ex');
			%x and y are as defined in the sample space. We need to generate x
			%and y in the detector space which may be done as

			if numel(x)>1
				x = -1*(x(1) - (0:(numel(x)-1))*diff(x(1:2)) + numel(x)*diff(x(1:2)));
			end

			if numel(y)>1
				y = -1*(y(1) - (0:(numel(y)-1))*diff(y(1:2)) + numel(y)*diff(y(1:2)));
			end

			[X,Y,Z] = ndgrid(x,y,z);

			%%%%%%%%%%%%%%%
			%xx = x;
			%yx = y;
			%%%%%%%%%%%%%%%

			eval(sprintf('Dx = %s(X,Y,beta_det,detmodevec(im));', detsensefun));

			if im==1
				Dx_tilde = zeros(numel(detmodevec),size(Dx,1),size(Dx,2));
			end
			Dx_tilde(im,[1 fliplr(2:end)],[1 fliplr(2:end)]) = ifft2(Dx);


			[x,y,z] = yeeposition(i_obs,j_obs,k_det_obs-illorigin(3)+Dzl,delta,'Ey');
			%x and y are as defined in the sample space. We need to generate x
			%and y in the detector space which may be done as

			if numel(x)>1
				x = -1*(x(1) - (0:(numel(x)-1))*diff(x(1:2)) + numel(x)*diff(x(1:2)));
			end

			if numel(y)>1
				y = -1*(y(1) - (0:(numel(y)-1))*diff(y(1:2)) + numel(y)*diff(y(1:2)));
			end

			z_obs = z;
			[X,Y,Z] = ndgrid(x,y,z);
			eval(sprintf('Dy = %s(X,Y,beta_det,detmodevec(im));', detsensefun));
			if im==1
				Dy_tilde = zeros(numel(detmodevec),size(Dy,1),size(Dy,2));
			end
			Dy_tilde(im,[1 fliplr(2:end)],[1 fliplr(2:end)]) = ifft2(Dy);
		end

		D_tilde.Dx_tilde = reshape(Dx_tilde,numel(detmodevec),size(Dy,1),size(Dy,2));
		D_tilde.Dy_tilde = reshape(Dy_tilde,numel(detmodevec),size(Dy,1),size(Dy,2));
		%save detvars Dx Dx_tilde fx_vec fy_vec xx yx i_obs j_obs;
    else
		Pupil = [];

		k_det_obs_global = [];

		f_vec.fx_vec = [];
		f_vec.fy_vec = [];

		D_tilde.Dx_tilde = [];
		D_tilde.Dy_tilde = [];
    end
    %finished seting up detector integral data, if required

    %save detdata Dy_tilde Dx_tilde Pupil fx_vec fy_vec;

    %begin iterations
    tic; %start the clock ticking
	 %for tind = 1:Nt

	tind_start = 1;
	tind_end   = 1; %the case of runmode = 'complete'
	if strncmp(runmode,'analyse',7)
	    tind_end = Nt;
	end

	if strncmp(operation,'run',3)
	    %first get the path
	    wd = pwd;
	    unify_match = strfind(wd,'dispersive');
	    dash_match = strfind(wd,'\');
	    if isempty(dash_match)
			dash_match = strfind(wd,'/');
	    end
	    [min_match,inds] = find(unify_match<dash_match);
	    addpath(sprintf('%siterater',wd(1:dash_match(inds(1)))));
	    fprintf(1,'Starting iterations\n');
	    tic;
	    for tind = tind_start:tind_end
			if(toc > 1)
				fprintf(1,'Done %d of %d\n',tind-tind_start,length(tind_start:tind_end));
				tic;
			end

			%Must ensure that the iwave_l values are complex
			if interface.I0(2) | interface.I1(2)
				Isource = complex(real(Isource),imag(Isource));
			end
			if interface.J0(2) | interface.J1(2)
				Jsource = complex(real(Jsource),imag(Jsource));
			end
			if interface.K0(2) | interface.K1(2)
				Ksource = complex(real(Ksource),imag(Ksource));
			end

			if strncmp(runmode,'analyse',7)
				%[Ex_out, Ey_out, Ez_out, Ex_bs, Ey_bs, Hx_bs, Hy_bs] =
				%iterater(Cmaterial,Dmaterial,C,D,freespace,interface,Isource,Jsource,Ksource,tvec_E,
				%tvec_H,omega_an,to_l,hwhm,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu,double(strcmp(lower_boundary_update,'true')),tind,dt,tind-1,sourcemode,runmode);
				[Ex_out, Ey_out, Ez_out, Hx_out, Hy_out, Hz_out, Ex_bs,Ey_bs, Hx_bs, Hy_bs, x_out, y_out, z_out,Ex_i, Ey_i, Ez_i, Hx_i, Hy_i, Hz_i,x_i,y_i,z_i,vertices,camplitudes,facets,maxresfield] = iterater(fdtdgrid,Cmaterial,Dmaterial,C,D,freespace,disp_params,delta,interface,Isource,Jsource,Ksource,grid_labels,omega_an,to_l,hwhm,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu,tind,dt,tind-1,sourcemode,runmode,exphasorsvolume,exphasorssurface,phasorsurface,phasorinc,dimension,conductive_aux,dispersive_aux,structure,f_ex_vec);
			else
				%[Ex_out, Ey_out, Ez_out, Ex_bs, Ey_bs, Hx_bs, Hy_bs] =
				%iterater(Cmaterial,Dmaterial,C,D,freespace,interface,Isource,Jsource,Ksource,tvec_E,
				%tvec_H,omega_an,to_l,hwhm,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu,double(strcmp(lower_boundary_update,'true')),Nt,dt,0,sourcemode,runmode);
				[Ex_out, Ey_out, Ez_out,  Hx_out, Hy_out, Hz_out, Ex_bs, Ey_bs, Hx_bs, Hy_bs, x_out, y_out, z_out,Ex_i, Ey_i, Ez_i, Hx_i, Hy_i, Hz_i,x_i,y_i,z_i,vertices,camplitudes,facets,maxresfield] = iterater(fdtdgrid,Cmaterial,Dmaterial,C,D,freespace,disp_params,delta,interface,Isource,Jsource,Ksource,grid_labels,omega_an,to_l,hwhm,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu,Nt,dt,0,sourcemode,runmode,exphasorsvolume,exphasorssurface,phasorsurface,phasorinc,dimension,conductive_aux,dispersive_aux,structure,f_ex_vec);
			end

			if strncmp(runmode,'analyse',7)
				for lvar=1:length(outputs_array)
					%t_output_matrix = 0;
					start_cell = outputs_array{lvar}{3};
					end_cell = outputs_array{lvar}{4};
					if ~isempty(findstr(outputs_array{lvar}{5},'x'))
						t_output_matrix.Ex = squeeze(fdtdgrid.Exy(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)) + fdtdgrid.Exz(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)));
						t_output_matrix.Hx = ...
							squeeze(fdtdgrid.Hxy(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)) + fdtdgrid.Hxz(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)));
					end
					if ~isempty(findstr(outputs_array{lvar}{5},'y'))
						t_output_matrix.Ey = squeeze(fdtdgrid.Eyx(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)) + fdtdgrid.Eyz(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)));
						t_output_matrix.Hy = squeeze(fdtdgrid.Hyx(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)) + fdtdgrid.Hyz(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)));
					end
					if ~isempty(findstr(outputs_array{lvar}{5},'z'))
						t_output_matrix.Ez = squeeze(fdtdgrid.Ezx(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)) + fdtdgrid.Ezy(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)));
						t_output_matrix.Hz = squeeze(fdtdgrid.Hzx(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)) + fdtdgrid.Hzy(start_cell(1):end_cell(1),start_cell(2):end_cell(2),start_cell(3):end_cell(3)));
					end
					if strcmp(outputs_array{lvar}{6},'dump')
						fr = t_output_matrix;
						fid = fopen(sprintf('%s.mat',outputs_array{lvar}{1}),'a');
						if fid ~= -1
							fclose(fid);
							save(sprintf('%s%d',outputs_array{lvar}{1},tind),'fr');
						else
							fprintf('Failed to open %s.mat for reading - skipping\n',outputs_array{lvar}{1});
						end
					else
						accumulate_output{lvar}{tind}= t_output_matrix;
					end
				end
			end
		end %for tind

		%save the accumulated field files
		if strncmp(runmode,'analyse',7)
			for lvar=1:length(outputs_array)
				if strcmp(outputs_array{lvar}{6},'accumulate')
					fr = accumulate_output{lvar};
					fid = fopen(sprintf('%s.mat',outputs_array{lvar}{1}),'a');
					if fid ~= -1
						fclose(fid);
						save(outputs_array{lvar}{1},'fr');
					else
						fprintf('Failed to open %s.mat for reading - skipping\n',outputs_array{lvar}{1});
					end
				end
			end
		end

		fprintf(1,'%f per iteration\n',toc/Nt);
	else
		fprintf(1,'Writing %s...',outfile);
		%we don't want to save the whole grid here just the materials
		materials = fdtdgrid.materials;
		clear fdtdgrid;
		fdtdgrid.materials = materials;
		%add in the information for allocating the field arrays
		fdtdgrid.I_tot = I_tot;
		fdtdgrid.J_tot = J_tot;
		fdtdgrid.K_tot = K_tot;

		%now save the file
		%iterater(fdtdgrid,Cmaterial,Dmaterial,C,D,freespace,interface,Isource,Jsource,Ksource,tvec_E,
		%tvec_H,omega_an,to_l,hwhm,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu,double(strcmp(lower_boundary_update,'true')),Nt,dt,0,sourcemode,runmode);
		%added the one to avoid a bug in casting to double -
			%clashed with release 13.

		tind = 0;

		if interface.I0(2) | interface.I1(2)
			Isource = complex(real(Isource),imag(Isource));
		end
		if interface.J0(2) | interface.J1(2)
			Jsource = complex(real(Jsource),imag(Jsource));
		end
		if interface.K0(2) | interface.K1(2)
			Ksource = complex(real(Ksource),imag(Ksource));
		end

		vers = version;
		campssample2 = campssample;
		if strncmp(operation,'filesetup',9)
			save(outfile,'fdtdgrid','Cmaterial','Dmaterial','C','D','freespace','interface','Isource','Jsource','Ksource','grid_labels','omega_an','to_l','hwhm','Dxl','Dxu','Dyl','Dyu','Dzl','Dzu','Nt','dt','tind','sourcemode','runmode','exphasorsvolume','exphasorssurface','intphasorssurface','phasorsurface','phasorinc','disp_params','delta','dimension','conductive_aux','dispersive_aux','structure','f_ex_vec','exdetintegral','f_vec','Pupil','D_tilde','k_det_obs_global','air_interface','intmatprops','use_bli','tdfield','tdfdir','fieldsample','campssample','use_pstd','-v7.3');
		else
			save(outfile,'fdtdgrid');
		end

		fprintf(1,'Done\n');
		%just to avoid a warning about setting output variables
		Ex_out = [];
		Ey_out = [];
		Ez_out = [];

		Hx_out = [];
		Hy_out = [];
		Hz_out = [];

		Ex_bs = [];
		Ey_bs = [];


		Hx_bs = [];
		Hy_bs = [];

		x_out = [];
		y_out = [];
		z_out = [];

		Ex_i = [];
		Ey_i = [];
		Ez_i = [];
		Hx_i = [];
		Hy_i = [];
		Hz_i = [];
		x_i = [];
		y_i = [];
		z_i = [];

		vertices = [];
		camplitudes = [];
		facets = [];
		maxresfield = [];
	end
end
end % function iteratefdtd_matrix

function result = has_ijk_source_matricies(data)
    fields = fieldnames(data);
    if ~(numel(fields) == 3)
        result = false;
    else
        result = strcmp(fields(1), 'Isource') & ...
                 strcmp(fields(2), 'Jsource') & ...
                 strcmp(fields(3), 'Ksource');
    end
end


function result = has_exi_eyi(data)
    fields = fieldnames(data);
    if ~(numel(fields) == 2)
        result = false;
    else
        result = strcmp(fields(1), 'exi') & ...
                 strcmp(fields(2), 'eyi');
    end
end


function assert_are_not_defined(efname, hfname)
    assert(strlength(efname) == 0, ...
        'TDMSException:IncompatibleInput', ...
        'An efield should not be defined. Set efname to an empty string');
    assert(strlength(hfname) == 0, ...
        'TDMSException:IncompatibleInput', ...
        'A hfield should not be defined. Set hfname to an empty string');
end

function assert_are_defined(efname, hfname)
    assert(strlength(efname) > 0, ...
        'TDMSException:IncompatibleInput', 'An efname must be defined');
    assert(strlength(hfname) > 0, ...
        'TDMSException:IncompatibleInput', 'A hfname must be defined');
end

function assert_source_has_correct_dimensions(Isource, Jsource, Ksource, interface)

    [mI,nI,oI] = size(Isource);
    [mJ,nJ,oJ] = size(Jsource);
    [mK,nK,oK] = size(Ksource);

    %Now make sure that the source matrices have the correct dimensions
    if length(Isource) > 0  % Check that the Isource dimensions are correct
        if ~( (mI==8) & (nI==(interface.J1(1) - interface.J0(1) + 1)) & (oI==(interface.K1(1) - interface.K0(1) + 1)))
            error('TDMSException:InvalidIlluminationDimensions',...
                  'Isource read in from %s might has incorrect dimenions');
        end
    end

    if length(Jsource) > 0  % Check that the Jsource dimensions are correct
        if ~( (mJ==8) & (nJ==(interface.I1(1) - interface.I0(1) + 1)) & (oJ==(interface.K1(1) - interface.K0(1) + 1)))
            error('TDMSException:InvalidIlluminationDimensions',...
                  'Jsource read in from %s might has incorrect dimenions');
        end
    end

    if length(Ksource) > 0  % Check that the Ksource dimensions are correct
        if ~( (mK==8) & (nK==(interface.I1(1) - interface.I0(1) + 1)) & (oK==(interface.J1(1) - interface.J0(1) + 1)))
            error('TDMSException:InvalidIlluminationDimensions',...
                  'Ksource read in from %s has incorrect dimenions');
        end
    end
end

function assert_exi_eyi_have_correct_dimensions(exi, eyi, I_tot, J_tot, Nt)

    [mX,nX,oX] = size(exi);
    [mY,nY,oY] = size(eyi);

    if ~(mX == I_tot + 1 & nX == J_tot + 1 & oX == Nt)
        error('TDMSException:InvalidIlluminationDimensions',...
              sprintf('exi must have dimensions (%d, %d, %d)', I_tot + 1, J_tot + 1, Nt));
    end

    if ~(mY == I_tot + 1 & nY == J_tot + 1 & oY == Nt)
        error('TDMSException:InvalidIlluminationDimensions',...
              sprintf('eyi must have dimensions (%d, %d, %d)', I_tot + 1, J_tot + 1, Nt));
    end
end

function interface_component_valid(C0, C1, upper, lower)
	%% Determine if the interface component in the C-plane is valid. Throw error if not.
	%% C represents either I, J, or K.
	% An interface component is valid when either:
	% 	C0(2) and C1(2) are both non zero AND:
	%		lower <= C0(1) <= C1(1) <= upper
	%	C0(2) OR C1(2) is zero
	if ~exists('lower', 'var')
		lower = 1;
	end

	if C0(2) && C1(2)
		if C0(1) < lower
			error('Interface component less than %d', lower);
		elseif C0(1) > C1(1)
			error('C0(1) must be less than C1(1)');
		elseif C1(1) > upper
			error('Interface component greater than %d', upper);
		end
	end
end
