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

%input the configuration information
[fid_input,message] = fopen(input_file,'r');

%check if file was opened successfully
if fid_input== -1
    error(sprintf('File %s could not be opened for reading',input_file));
end

%proceed to_l read in config information
current_line = fgets(fid_input);
while current_line ~= -1
    if ~isempty(findstr('material_file',current_line))
	if isempty(material_file) | ((~isempty(material_file)) & (findstr('material_file',current_line)>findstr('%',current_line)))
	    eval(current_line);
	else
	    fprintf(1,'material_file specified as argument and in %s, using %s\n',input_file,material_file);
	end
    else
	eval(current_line);
    end

    current_line = fgets(fid_input);
end

if isempty(material_file)
    clear material_file;
end
%now need to_l check that all of the required variables have been set
variables = {'delta','I','J','K','n','R0','Dxl','Dxu','Dyl','Dyu','Dzl','Dzu','dt','epsr','mur','f_an','Nt','interface','material_file','efname','hfname','wavelengthwidth','z_launch','illorigin','runmode','sourcemode','exphasorsvolume','exphasorssurface','intphasorssurface','phasorsurface','phasorinc','dimension','multilayer','kappa_max','vc_vec','wp_vec','structure','f_ex_vec','exdetintegral','k_det_obs','NA_det','beta_det','detmodevec','detsensefun','air_interface','intmatprops','intmethod','tdfdir','fieldsample','campssample','usecd','compactsource'};
must_abort = 0; %assumes all variables have been defined
for lvar = 1:length(variables)
    if exist(variables{lvar}) ~= 1
	if strncmp(variables(lvar),'z_launch',8)
	    fprintf(1,'Failed to define %s, setting it to 0\n',variables{lvar});
	    z_launch = 0;
	elseif strncmp(variables(lvar),'runmode',7)
	    fprintf(1,'Failed to define %s, setting it to ''complete''\n',variables{lvar});
	    runmode = 'complete';
	elseif strncmp(variables(lvar),'dimension',9)
	    fprintf(1,'Failed to define %s, setting it to ''3''\n',variables{lvar});
	    dimension = '3';
	elseif strncmp(variables(lvar),'phasorinc',9)
	    fprintf(1,'Failed to define %s, setting it to [1 1 1]\n',variables{lvar});
	    phasorinc = [1 1 1];
	elseif strncmp(variables(lvar),'multilayer',10)
	    fprintf(1,'Failed to define %s, setting it to []\n',variables{lvar});
	    multilayer = [];
	elseif strncmp(variables(lvar),'kappa_max',9)
	    fprintf(1,'Failed to define %s, setting it to 1\n',variables{lvar});
	    kappa_max = 1;
	elseif strncmp(variables{lvar},'vc_vec',6)
	    tmp_str = num2str(zeros(1,length(epsr)));
	    fprintf(1,'Failed to define %s, setting it to [%s]\n',variables{lvar},tmp_str);
	    vc_vec = zeros(1,length(epsr));
	elseif strncmp(variables{lvar},'wp_vec',6)
	    tmp_str = num2str(zeros(1,length(epsr)));
	    fprintf(1,'Failed to define %s, setting it to [%s]\n',variables{lvar},tmp_str);
	    wp_vec = zeros(1,length(epsr));
	elseif strncmp(variables{lvar},'structure',9)
	    fprintf(1,'Failed to define %s, setting it to []\n',variables{lvar});
	    structure = [];
	elseif strncmp(variables{lvar},'f_ex_vec',8)
	    fprintf(1,'Failed to define %s, setting it to f_an\n',variables{lvar});
	    f_ex_vec = f_an;
	elseif strncmp(variables{lvar},'intphasorssurface',17)
	    fprintf(1,'Failed to define %s, setting it to 1\n',variables{lvar});
	    intphasorssurface=1;
	elseif strncmp(variables{lvar},'exdetintegral',13)
	    fprintf(1,'Failed to define %s, setting it to 0\n',variables{lvar});
	    exdetintegral = 0;
	elseif strncmp(variables{lvar},'k_det_obs',9)
	    if exdetintegral == 1
		fprintf(1,'Failed to define %s, yet exdetintegral is set to 1\n',variables{lvar});
		must_abort = 1;
	    end
	elseif strncmp(variables{lvar},'NA_det',6)
	    if exdetintegral == 1
		fprintf(1,'Failed to define %s, yet exdetintegral is set to 1\n',variables{lvar});
		must_abort = 1;
	    end
	elseif strncmp(variables{lvar},'beta_det',8)
	    if exdetintegral == 1
		fprintf(1,'Failed to define %s, yet exdetintegral is set to 1\n',variables{lvar});
		must_abort = 1;
	    end
	elseif strncmp(variables{lvar},'detmodevec',10)
	    if exdetintegral == 1
		fprintf(1,'Failed to define %s, yet exdetintegral is set to 1\n',variables{lvar});
		must_abort = 1;
	    end
	elseif strncmp(variables{lvar},'detsensefun',11)
	    if exdetintegral == 1
		fprintf(1,'Failed to define %s, yet exdetintegral is set to 1\n',variables{lvar});
		must_abort = 1;
	    end
	elseif strncmp(variables{lvar},'air_interface',13)
	    fprintf(1,'Failed to define %s, setting it to []\n',variables{lvar});
	    air_interface = [];
	elseif strncmp(variables{lvar},'intmatprops',13)
	    fprintf(1,'Failed to define %s, setting it to 1\n',variables{lvar});
	    intmatprops = 1;
	elseif strncmp(variables{lvar},'intmethod',9)
	    fprintf(1,'Failed to define %s, setting it to 1\n',variables{lvar});
	    intmethod = 1;
	elseif strncmp(variables{lvar},'tdfdir',6)
	    fprintf(1,'Failed to define %s, setting it to empty string\n',variables{lvar});
	    tdfdir = '';
	elseif strncmp(variables{lvar},'fieldsample',11)
	    fprintf(1,'Failed to define %s, setting it to []\n',variables{lvar});
	    fieldsample.i = [];
	    fieldsample.j = [];
	    fieldsample.k = [];
	    fieldsample.n = [];
	elseif strncmp(variables{lvar},'campssample',11)
	    fprintf(1,'Failed to define %s, setting campssample.vertices = [] and campssample.components = []\n',variables{lvar});
	    campssample.vertices = [];
	    campssample.components = [];
	elseif strncmp(variables{lvar},'usecd',5)
	    fprintf(1,'Failed to define %s, setting it to 1 (meaning the FDTD algorithm will be used)\n',variables{lvar});
	    usecd=1;
	elseif strncmp(variables{lvar},'compactsource',13)
	    fprintf(1,'Failed to define %s, setting it to 0 (meaning that both efname and hfname should be specified)',variables{lvar});
	    compactsource=0;
	else
	    fprintf(1,'Failed to define %s\n',variables{lvar});
	    must_abort = 1;
	end
    end
end

if must_abort
    error('Not all variables were defined');
end

%just check that the outputs_array is valid
for lvar = 1:size(outputs_array,1)
    if length(outputs_array{lvar}) < 6
        error(sprintf('insufficient number of entries in outputs_array number %d',lvar));
    end
    if strcmp(outputs_array{lvar}{2},'interior') %means indexing is relative to the first non-pml cell
						 %must check that all nodes are within the allowed range
    node1 = outputs_array{lvar}{3};
    node2 = outputs_array{lvar}{4};
    if ~(node1(1)<=node2(1) & node1(2)<=node2(2) & node1(3)<=node2(3))
	error(sprintf('outputs_array entry %d range vectors incorrect',lvar));
    end

    if strcmp('3',dimension)
        if sum(node1<=[I J K] & node1>=[1 1 1] & node2<=[I J K] & node2>=[1 1 1])~=3
            error(sprintf('outputs_array entry %d range vectors incorrect',lvar));
        end
    else
        if sum(node1<=[I J 1] & node1>=[1 1 1] & node2<=[I J 1] & node2>=[1 1 1])~=3
            error(sprintf('outputs_array entry %d range vectors incorrect',lvar));
        end
    end

    %now set the values for global indexing
    outputs_array{lvar}{3} = outputs_array{lvar}{3} + [Dxl Dyl Dzl];
    outputs_array{lvar}{4} = outputs_array{lvar}{4} + [Dxl Dyl Dzl];

    elseif strcmp(outputs_array{lvar}{2},'global') %means indexing is relative to the first cell in the grid
        node1 = outputs_array{lvar}{3};
        node2 = outputs_array{lvar}{4};
        if ~(node1(1)<=node2(1) & node1(2)<=node2(2) & node1(3)<=node2(3))
            error(sprintf('outputs_array entry %d range vectors incorrect',lvar));
        end

        if sum(node1<=[(I+Dxl+Dxu+1) (J+Dyl+Dyu+1) (K+Dyl+Dyu+1)] & node1>=[1 1 1] & node2<=[(I+Dxl+Dxu+1) (J+Dyl+Dyu+1) (K+Dyl+Dyu+1)] & node2>=[1 1 1])~=3
            error(sprintf('outputs_array entry %d range vectors incorrect',lvar));
        end
    else
        error(sprintf('Indexing for outputs_array entry number %d is incorrect, should be interior or global',lvar));
    end
end

%check campssample is valid
if ~isempty(campssample)
    campssample_fieldnames = fieldnames(campssample);
    if numel( campssample_fieldnames ) ~= 2
	error(sprintf('campssample has been incorrectly specified\n'));
    end
    if ~((strcmp(campssample_fieldnames(2),'components') & strcmp(campssample_fieldnames(1),'vertices')) | (strcmp(campssample_fieldnames(1),'components') & strcmp(campssample_fieldnames(2),'vertices')))
	error(sprintf('campssample has been incorrectly specified\n'));
    end
    campssample.vertices = int32( campssample.vertices );
    campssample.components = int32( campssample.components );
end

%just check that the outputs_array is valid - in particular the components
for lvar = 1:length(outputs_array)
    if length(outputs_array{lvar}{5}) > 0
        for mvar=1:length(outputs_array{lvar}{5})
            if isempty(findstr(outputs_array{lvar}{5}(mvar),'xyz'))
                error(sprintf('%s is not a vector component',  outputs_array{lvar}{5}(mvar)));
            end
        end
    else
        error(sprintf('no components were specified in outputs_array element %d',lvar));
    end
end

%just check that the outputs_array is valid - in particular that the output style is accumulate or dump
for lvar = 1:length(outputs_array)
    if ~( strcmp(outputs_array{lvar}{6},'accumulate') | strcmp(outputs_array{lvar}{6},'dump') )
        error(sprintf('mode of output writing is incorrect for outputs_array entry %d',lvar));
    end
end

%check that interface has been fully specified, if not, any
%un-specified interface conditions will be set to 0, ie, the
%interface condition will not be implemented
interface_field = {'I0','I1','J0','J1','K0','K1'};
%global Isource Jsource Ksource interface source_field x y z X Y Z;
for lvar = 1:length(interface_field)
    if ~isfield(interface,interface_field{lvar})
	fprintf('interface.%s has not been defined, setting it to [0 0]\n', interface_field{lvar});
	interface = setfield(interface,interface_field{lvar},[0 0]);
    end
end

%now check that the entries for I0...K1 are valid
if interface.I0(2) & interface.I1(2)
    if interface.I1(1) < interface.I0(1)
	error('Must have interface.I1>=interface.I0');
    end
    if interface.I1(1) > I
	error('Must have interface.I1 <= I');
    end
    if interface.I0(1) < 1
	error('Must have interface.I0 >= 1');
    end
end

if interface.J0(2) & interface.J1(2)
    if interface.J1(1) < interface.J0(1)
	error('Must have interface.J1>=interface.J0');
    end
    if interface.J1(1) > J
	error('Must have interface.J1 <= J');
    end
    if interface.J0(1) < 1
	error('Must have interface.J0 >= 1');
    end
end

if interface.K0(2) & interface.K1(2)
    if interface.K1(1) < interface.K0(1)
	error('Must have interface.K1>=interface.K0');
    end
    if interface.K1(1) > K
	error('Must have interface.K1 <= K');
    end
    if interface.K0(1) < 1
	error('Must have interface.K0 >= 1');
    end
end

%check that runmode is valid
if ~(strncmp(runmode,'analyse',7) | strncmp(runmode,'complete',8))
    error(sprintf('runmode ''%s'' invalid',runmode));
end

%check that sourcemode is valid
if ~(strncmp(sourcemode,'pulsed',6) | strncmp(sourcemode,'steadystate',11))
    error(sprintf('sourcemode ''%s'' invalid',sourcemode));
end

%check that operation is valid
if ~(strncmp(operation,'run',3) | strncmp(operation,'filesetup',9) | strncmp(operation,'gridsetup',9) | strncmp(operation,'illsetup',8))
    error(sprintf('operation ''%s'' invalid',operation));
elseif strncmp(operation,'filesetup',9)
    if strncmp(runmode,'analyse',7)
	error('Only runmode complete is suported for operation filesetup');
    end
end

%check that phasorsurface is valid
if exphasorssurface
    [psm, psn] = size(phasorsurface);
    if psm*psn ~= 6 %check for correct number of elements
	error('phasorsurface must be a vector of 6 elements');
    else            %now check that each entry is valid
	if phasorsurface(1) + Dxl -2 < 1
	    error('phasorsurface(1) out of range');
	elseif phasorsurface(2) - Dxu > I
	    error('phasorsurface(2) out of range');
	elseif (phasorsurface(3) + Dyl -2 < 1) & ((J + Dyl + Dyu)>0)
	    error('phasorsurface(3) out of range');
	elseif (phasorsurface(4) - Dyu > J) & ((J + Dyl + Dyu)>0)
	    error('phasorsurface(4) out of range');
	elseif (phasorsurface(5) + Dzl -2 < 1) & ~( (K==0) & (phasorsurface(5)==0) )
	    error('phasorsurface(5) out of range');
	elseif phasorsurface(6) - Dzu > K
	    error('phasorsurface(6) out of range');
	end


% $$$     if phasorsurface(1) < 1
% $$$ 	error('phasorsurface(1) out of range');
% $$$     elseif phasorsurface(2) > I
% $$$ 	error('phasorsurface(2) out of range');
% $$$     elseif phasorsurface(3) < 1
% $$$ 	error('phasorsurface(3) out of range');
% $$$     elseif phasorsurface(4) > J
% $$$ 	error('phasorsurface(4) out of range');
% $$$     elseif (phasorsurface(5) < 1) & ~( (K==0) & (phasorsurface(5)==0) )
% $$$ 	error('phasorsurface(5) out of range');
% $$$     elseif phasorsurface(6) > K
% $$$ 	error('phasorsurface(6) out of range');
% $$$     end
% $$$     %now convert to global coordinates
	phasorsurface = phasorsurface + [Dxl Dxl Dyl Dyl Dzl Dzl];
	%still in the matlab convention though

    end
end

%check to ensure that the phasor surface has been specified
%correctly
if ( mod( (phasorsurface(2) - phasorsurface(1)),phasorinc(1) ) | mod( (phasorsurface(4) - phasorsurface(3)),phasorinc(2) ) | mod( (phasorsurface(6) - phasorsurface(5)),phasorinc(3) ) )
    mod( (phasorsurface(2) - phasorsurface(1)),phasorinc(1) )
    mod( (phasorsurface(4) - phasorsurface(3)),phasorinc(2) )
    mod( (phasorsurface(6) - phasorsurface(5)),phasorinc(3) )
    error('incorrect specification of phasorinc');
end

%check that intmethod has valid values, ie, either 1 or 2
%1 = cubic
%2 = band limited
if ~( intmethod==1 | intmethod==2)
    error('incorrect specification of intmethod');
end


%this is really just for completeness
x_max = I*delta.x;
y_max = J*delta.y;
z_max = K*delta.z;

%now check that dimension is correct
if ~(strcmp('3',dimension) | strcmp('TE',dimension) | strcmp('TM',dimension))
    error('Dimension must take the value ''3'', ''TM'' or ''TE''');
end

%check that multilayer and epsr have the correct dimension
if isempty(multilayer)
    if numel(epsr)~=1
	error('epsr should have only a single element in the case of a single layer');
    end
else
    if numel(multilayer)~=(numel(epsr)-1)
	error('epsr should have one more member than multilayer');
    end
end

I_tot = I + Dxl + Dxu;  %Extent of the grid - including the PML
J_tot = J + Dyl + Dyu;
K_tot = K + Dzl + Dzu;

%convert campssample.vertices to global coordinates, matlab version
%not C
if ~isempty(campssample.vertices)
    campssample.vertices(:,1) = campssample.vertices(:,1) + Dxl;
    campssample.vertices(:,2) = campssample.vertices(:,2) + Dyl;
    campssample.vertices(:,3) = campssample.vertices(:,3) + Dzl;
else
    campssample.components = [];%just make sure
end

%error checking for structure
if isempty(multilayer)
    if ~isempty(structure)
	error('structure is not empty yet multlayer is');
    end
else
    if  ~isempty(structure)
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
	if  ~tmp_pro(end);
	    error('structure should have first row monotonically increasing');
	end

	if structure(1,1)<1
	    error('structure should have first row entries all greater than 1');
	end

	if structure(1,end)>(I_tot+1)
	    error('structure should have first row entries all less than I+1');
	end

	[ml_mat,pl_mat] = ndgrid(multilayer,structure(2,:));
	if ~isempty(find( (ml_mat+pl_mat)<2 ))
	    error('trench is breaking PML boundary, reduce magnitude of structure displacement');
	end

	if ~isempty(find( (ml_mat+pl_mat)>(K-1) ))
	    error('trench is breaking PML boundary, reduce magnitude of structure displacement');
	end
	%%
	i_int = 1:(I_tot+1);
	p_int = int32(round(interp1(structure(1,:),structure(2,:),i_int,'linear')));
	i_int = int32(i_int);

	structure = [i_int;p_int];
	%%
    end
end
%error checking for structure

%check that kappa_max has the correct dimension
if numel(kappa_max)~=1
    if numel(kappa_max)~=numel(epsr)
	error('should have numel(kappa_ax)=numel(epsr)');
    end
end

%check that vc_vec and wp_vec have the correct dimension
if numel(vc_vec) ~= numel(wp_vec)
    error('vc_vec and wp_vec should have the same number of elements');
end

if numel(multilayer)~=(numel(vc_vec)-1)
    error('vc_vec and wp_vec should have one more member than multilayer');
end

%check correctness of multilayer
if ~isempty(multilayer)
    if min(multilayer) <1
	error('multilayer must have elements greater than 0');
    end
    if max(multilayer)>K
	error('multilayer must have elements less than or equal to  K');
    end
    if sum(sort(multilayer)==multilayer)~=numel(multilayer)
	error('multilayer must be monotonically increasing');
    end
end

%check the correctness of f_ex_vec
if numel(f_ex_vec)>1
    if ~(size(f_ex_vec,1)==1 | size(f_ex_vec,2)==1)
	error('f_ex_vec should be a vector (ie, a matrix with one singleton dimension)');
    end
end

%check correctness of detector sensitivity evaluation parameters
if ~(exdetintegral==0 | exdetintegral==1)
    error('exdetintegral should take the value 0 or 1');
end

if exdetintegral==1
    if NA_det<0 | NA_det>1
	error('NA_det should be between 0 and 1');
    end

    if sum(abs( round(detmodevec) - detmodevec ))~=0
	error('detmodevec should contain only integers');
    end

    if ~all( detmodevec > 0)
	error('detmodevec should contain only positive integers');
    end

    if numel(detmodevec)==0
	error('detmodevec should contain at least one element');
    end
    detmodevec = detmodevec(:);
end

%check that tdfdir does not end with a /
if(~isempty(tdfdir))
    if tdfdir(end)=='/'
	tdfdir = tdfdir(1:(end-1));
    end
end

%check that air_interface is correct
if ~isempty(air_interface)
    if air_interface > -illorigin(3)*delta.z
	error(sprintf('\nError: air_interface is intended to be defined\noutside of the simulation space (z=%e),\nin the negative z-direction, whereas it is at\nz=%e\n',-illorigin(3)*delta.z,air_interface));
    end
end

%here we need to establish if exi and eyi have been passed
exipresent=0;
if length(ill_file) > 0%must have already computed the illumination source
    data = load(ill_file);
    %here we can have a data file with elemenets Isource, Jsource
    %and Ksource *or* exi and eyi
    fieldnames_ill = fieldnames(data);
    clear data;
    if numel(fieldnames_ill)==2
	exipresent=1;
    end
end
fprintf(1,'exipresent: %d\n',exipresent);

%check that compactsource, usecd and exipresent have compatible
%values

if ~usecd & ~compactsource
    error('A compact source must be used when not using central differences');
elseif compactsource & ~isempty(hfname)
    error('When using a compact source, the magnetic field input file must be empty');
elseif ~exipresent & compactsource & isempty(efname)
    error('When using a compact source, the electric field input file must be non-empty');
elseif ~exipresent & usecd & ~compactsource & ~(~isempty(efname) & ~isempty(hfname))
    error('When using central differences and a non-compact source condition, both the electric and magnetic field input files must be non-empty');
elseif exipresent & usecd & ~compactsource & (~isempty(efname) & isempty(hfname) | isempty(efname) & ~isempty(hfname))
    error('When specifying exi and/or eyi, along with a non-compact source, the electric and magnetic field input files should both be non-empty or both empty');
end

fprintf('Allocating grid...');
fdtdgrid = initialisesplitgrid(I,J,K,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu);
fprintf('Done\n');

[epso , muo , c] = import_constants;

%We now implement an incident planewave using a total/scattered field
%formulation. For k<K1 we have scattered and for k>=K1 we have total.
%The field is coupled into system by means of update equations to_l maintain
%consistency.

%correct:
refractive_index = sqrt(real(epsr(1)));

omega_an = 2*pi*f_an;
lambda_an = c/(f_an*refractive_index);
wave_num_an = 2*pi/lambda_an;%wave number in m^-1

fprintf('Initialising source field...');

%Isource(:,1,1) = [Ey Ez Hy Hz Ey Ez Hy Hz]
%Jsource(:,1,1) = [Ex Ez Hx Hz Ex Ez Hx Hz]
%Ksource(:,1,1) = [Ex Ey Hx Hy Ex Ey Hx Hy]

%now need to adjust interface to be in global coordinates
interface.I0(1) = interface.I0(1) + Dxl;
interface.I1(1) = interface.I1(1) + Dxl;
interface.J0(1) = interface.J0(1) + Dyl;
interface.J1(1) = interface.J1(1) + Dyl;
interface.K0(1) = interface.K0(1) + Dzl;
interface.K1(1) = interface.K1(1) + Dzl;

illorigin = illorigin + [Dxl Dyl Dzl];

%however, if the sourcemode is pulsed then we must reset these
%values to:


if strncmp(sourcemode,'pulsed',6)
    interface.I0(1) = 1;
    interface.J0(1) = 1;
    interface.I1(1) = I_tot+1;
    interface.J1(1) = J_tot+1;

    interface.I0(2) = 0;
    interface.I1(2) = 0;
    interface.J0(2) = 0;
    interface.J1(2) = 0;
    interface.K1(2) = 0;

    if (interface.K0(2)==0) & K~=0
	error('Running in pulsed mode with k0[0]=0, there is no point running');
    end
end

if length(ill_file) > 0%must have already computed the illumination source
    data = load(ill_file);
    %here we can have a data file with elemenets Isource, Jsource
    %and Ksource *or* exi and eyi
    fieldnames_ill = fieldnames(data);
    if numel(fieldnames_ill)==3%Case Isource, Jsource and Ksource specified
	Isource = data.Isource;
	Jsource = data.Jsource;
	Ksource = data.Ksource;
	[mI,nI,oI] = size(Isource);
	[mJ,nJ,oJ] = size(Jsource);
	[mK,nK,oK] = size(Ksource);
	tdfield.exi = [];
	tdfield.eyi = [];

	%Now make sure that the source matrices have the correct
	%dimensions
	if ~( (mI==8) & (mJ==8) & (mK==8) & (nI==(interface.J1(1) - interface.J0(1) + 1)) & (nJ==(interface.I1(1) - interface.I0(1) + 1)) & (nK==(interface.I1(1) - interface.I0(1) + 1)) & (oI==(interface.K1(1) - interface.K0(1) + 1)) & (oJ==(interface.K1(1) - interface.K0(1) + 1)) & (oK==(interface.J1(1) - interface.J0(1) + 1)))
	    (fprintf(1,'Illumination matrices read in from %s might have incorrect dimenions',ill_file));
	end
    elseif numel(fieldnames_ill)==2%Case exi and eyi specified
%	exi = data.exi;
%	eyi = data.eyi;
	tdfield = data;
	if (interface.I0(2) | interface.I1(2)) & (~isempty(efname))% hfname may be empty or non-empty  & (~isempty(hfname))
	    Isource = zeros(8,interface.J1(1) - interface.J0(1) + 1, interface.K1(1) - interface.K0(1) + 1);
	else
	    Isource=[];
	end

	if (interface.J0(2) | interface.J1(2)) & (~isempty(efname))% hfname may be empty or non-empty  & (~isempty(hfname))
	    Jsource = zeros(8,interface.I1(1) - interface.I0(1) + 1, interface.K1(1) - interface.K0(1) + 1);
	else
	    Jsource = [];
	end

	if (interface.K0(2) | interface.K1(2)) & (~isempty(efname))% hfname may be empty or non-empty  & (~isempty(hfname))
	    Ksource = zeros(8,interface.I1(1) - interface.I0(1) + 1, interface.J1(1) - interface.J0(1) + 1);
	else
	    Ksource = [];
	end
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%If the user has specified field function names, they can
        %also specify a pulsed field.
	%

	if (~isempty(efname))% hfname may be empty or non-empty  & (~isempty(hfname))
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

		if (~isempty(hfname))
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

		if (~isempty(hfname))
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

		if (~isempty(hfname))
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

		if (~isempty(hfname))
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
		if ~compactsource
		    eval(sprintf('source_field = %s(X,Y,Z);',efname ));
		    Ksource(1,:,:) = source_field{1};
		else
		    eval(sprintf('source_field = %s(X,Y,Z-delta.z/2);',efname ));
		    Ksource(1,:,:) = 2*source_field{1};
		end


		%Ey, K0
		[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
		z = z + z_launch;
		[X,Y,Z] = ndgrid(x,y,z);
		if ~compactsource
		    eval(sprintf('source_field = %s(X,Y,Z);',efname ));
		    Ksource(2,:,:) = source_field{2};
		else
		    eval(sprintf('source_field = %s(X,Y,Z-delta.z/2);',efname ));
		    Ksource(2,:,:) = 2*source_field{2};
		end


		if (~isempty(hfname))
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

		if (~isempty(hfname))
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
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
else%an illumination file has not been specified
    tdfield.exi = [];tdfield.eyi = [];
    if (interface.I0(2) | interface.I1(2))  & (~isempty(efname))% & (~isempty(hfname))
	Isource = zeros(8,interface.J1(1) - interface.J0(1) + 1, interface.K1(1) - interface.K0(1) + 1);
    else
	Isource=[];
    end

    if (interface.J0(2) | interface.J1(2))  & (~isempty(efname))% & (~isempty(hfname))
	Jsource = zeros(8,interface.I1(1) - interface.I0(1) + 1, interface.K1(1) - interface.K0(1) + 1);
    else
	Jsource = [];
    end

    if (interface.K0(2) | interface.K1(2))  & (~isempty(efname))% & (~isempty(hfname))
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


    %Ey, I0
    if interface.I0(2) & (~isempty(efname))% & (~isempty(hfname))
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

	if ~isempty(hfname)
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
    end

    i_source = interface.I1(1) - illorigin(1);
    %Ey, I1
    if interface.I1(2) & (~isempty(efname))% & (~isempty(hfname))
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

	if ~isempty(hfname)
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
    end

    %Set up the Jsource field. This has to be defined on a 2d array
    %over the range (I0,I1)x(K0,K1). We calculate the field values
    %assuming an origin for the illumination
    i_source = (interface.I0(1):interface.I1(1)) - illorigin(1);
    j_source = interface.J0(1) - illorigin(2);
    k_source = (interface.K0(1):interface.K1(1))- illorigin(3);


    %Ey, J0
    if interface.J0(2) & (~isempty(efname))% & (~isempty(hfname))
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

	if ~isempty(hfname)
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
    end

    j_source = interface.J1(1) - illorigin(2);
    %Ey, J1
    if interface.J1(2) & (~isempty(efname))% & (~isempty(hfname))
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

	if ~isempty(hfname)
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
    end

    %Set up the Ksource field. This has to be defined on a 2d array
    %over the range (I0,I1)x(J0,J1). We calculate the field values
    %assuming an origin for the illumination
    i_source = (interface.I0(1):interface.I1(1)) - illorigin(1);
    j_source = (interface.J0(1):interface.J1(1)) - illorigin(2);
    k_source = interface.K0(1) - illorigin(3);


    %Ex, K0
    if interface.K0(2) & (~isempty(efname))% & (~isempty(hfname))

	[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
	z = z + z_launch;
	%fprintf(1,'%d %d %d %e\n',i_source,j_source,k_source,z);
	[X,Y,Z] = ndgrid(x,y,z);
	if ~compactsource
	    eval(sprintf('source_field = %s(X,Y,Z);',efname ));
	    Ksource(1,:,:) = source_field{1};
	else
	    eval(sprintf('source_field = %s(X,Y,Z-delta.z/2);',efname ));
	    Ksource(1,:,:) = 2*source_field{1};
	end


	%Ey, K0
	[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
	z = z + z_launch;
	[X,Y,Z] = ndgrid(x,y,z);
	if ~compactsource
	    eval(sprintf('source_field = %s(X,Y,Z);',efname ));
	    Ksource(2,:,:) = source_field{2};
	else
	    eval(sprintf('source_field = %s(X,Y,Z-delta.z/2);',efname ));
	    Ksource(2,:,:) = 2*source_field{2};
	end


	if ~isempty(hfname)
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
    end
    k_source = interface.K1(1) - illorigin(3);
    %Ex, K1
    if interface.K1(2) & (~isempty(efname)) & (~isempty(hfname))
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

	if ~isempty(hfname)
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
fprintf('Done\n');

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
    if 0
	for i=interface.I0(1):interface.I1(1)
	    for j=interface.J0(1):interface.J1(1)
		for k=interface.K0(1):interface.K1(1)
		    i_source = i-illorigin(1);
		    j_source = j-illorigin(2);
		    k_source = k-illorigin(3);
		    if i<interface.I1(1)
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
			eval(sprintf('source_field = %s(x,y,z);',efname ));
			fdtdgrid.Exz(i,j,k) = real(source_field{1}*phasetermE);
		    end
		    if j<interface.J1(1)
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
			eval(sprintf('source_field = %s(x,y,z);',efname ));
			fdtdgrid.Eyz(i,j,k) = real(source_field{2}*phasetermE);
			if i==interface.I0(1)
			    if abs(source_field{2}-Isource(1,j-interface.J0(1)+1,k-interface.K0(1)+1))>1e-10
				fprintf(1,'Error - Isource\n');
			    end
			end
		    end
		    if k<interface.K1(1)
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ez');
			eval(sprintf('source_field = %s(x,y,z);',efname ));
			fdtdgrid.Ezx(i,j,k) = real(source_field{3}*phasetermE);
		    end
		    if (j<interface.J1(1)) & (k<interface.K1(1))
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hx');
			eval(sprintf('source_field = %s(x,y,z);',hfname ));
			fdtdgrid.Hxz(i,j,k) = real(source_field{1}*phasetermH);
		    end
		    if (i<interface.I1(1)) & (k<interface.K1(1))
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hy');
			eval(sprintf('source_field = %s(x,y,z);',hfname ));
			fdtdgrid.Hyz(i,j,k) = real(source_field{2}*phasetermH);
		    end
		    if (i<interface.I1(1)) & (j<interface.J1(1))
			[x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Hz');
			eval(sprintf('source_field = %s(x,y,z);',hfname ));
			fdtdgrid.Hzx(i,j,k) = real(source_field{3}*phasetermH);
		    end
		end
	    end
	end
    end

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
	 %for tind = 1:1

	 %for tind = 1:Nt

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
	     end %of tind?
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
		 rmpath('/home/ptpc2/prmunro/code/ptws1/FDTD/dispersive1.1/iterater');
		% rmpath('iterater');

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
		 save(outfile,'fdtdgrid','Cmaterial','Dmaterial','C','D','freespace','interface','Isource','Jsource','Ksource','grid_labels','omega_an','to_l','hwhm','Dxl','Dxu','Dyl','Dyu','Dzl','Dzu','Nt','dt','tind','sourcemode','runmode','exphasorsvolume','exphasorssurface','intphasorssurface','phasorsurface','phasorinc','disp_params','delta','dimension','conductive_aux','dispersive_aux','structure','f_ex_vec','exdetintegral','f_vec','Pupil','D_tilde','k_det_obs_global','air_interface','intmatprops','intmethod','tdfield','tdfdir','fieldsample','campssample','usecd','-v7.3');%,'-V6');%'-v7.3');%,'-V6');
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
