%function [ex_coords, ey_coords, tvec_E, fvec_E, f_an, hwhm, to_l] = getsourcecoords(input_file)
%
%input_file - file with input configuration information
%
function [ex_coords, ey_coords, tvec_E, fvec_E, f_an, hwhm, to_l] = getsourcecoords(input_file)

[fid_input,message] = fopen(input_file,'r');


%%controls whether to use FDTD or PSTD
useCD=0;
%%
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

material_file = [];

%now need to_l check that all of the required variables have been set
variables = {'delta','I','J','K','n','R0','Dxl','Dxu','Dyl','Dyu','Dzl','Dzu','dt','epsr','mur','f_an','Nt','interface','material_file','efname','hfname','wavelengthwidth','z_launch','illorigin','runmode','sourcemode','exphasorsvolume','exphasorssurface','phasorsurface','phasorinc','dimension','multilayer','kappa_max','vc_vec','wp_vec','structure','f_ex_vec','exdetintegral','k_det_obs','NA_det','detsensefun','beta_det','det_trans_x','det_trans_y','illspecfun'};
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
	elseif strncmp(variables{lvar},'exdetintegral',13)
	    fprintf(1,'Failed to define %s, setting it to 0\n',variables{lvar});
	    exdetintegral=0;
	elseif strncmp(variables(lvar),'k_det_obs',9)
	    fprintf(1,'Failed to define %s\n',variables{lvar});
	    must_abort=1;
	elseif strncmp(variables(lvar),'NA_det',6)
	    fprintf(1,'Failed to define %s\n',variables{lvar});
	    must_abort=1;
	elseif strncmp(variables(lvar),'beta_det',8)
	    fprintf(1,'Failed to define %s\n',variables{lvar});
	    must_abort=1;
	elseif strncmp(variables(lvar),'detsensefun',11)
	    fprintf(1,'Failed to define %s\n',variables{lvar});
	    must_abort=1;
	elseif strncmp(variables(lvar),'det_trans_x',11)
	    fprintf(1,'Failed to define %s, seeting to det_trans_x = 0\n',variables{lvar});
	    det_trans_x=0;
	elseif strncmp(variables(lvar),'det_trans_y',11)
	    fprintf(1,'Failed to define %s, setting to det_trans_y = 0\n',variables{lvar});
	    det_trans_y=0;
	elseif strncmp(variables(lvar),'illspecfun',10)
	    fprintf(1,'Failed to define %s, setting to illspecfun=''''\n',variables{lvar});
	    illspecfun = '';
	else
	    fprintf(1,'Failed to define %s\n',variables{lvar});
	    must_abort = 1;
	end
    end
end
k_obs = k_det_obs;
NA = NA_det;

%now check
% k_obs
% NA
% detsensefun

if exdetintegral==1
    for lvar = 1:length(variables)
	if exist(variables{lvar}) ~= 1
	    if strncmp(variables(lvar),'k_obs',4)
		fprintf(1,'Failed to define %s\n',variables{lvar});
		must_abort=1;
	    elseif strncmp(variables(lvar),'NA',2)
		fprintf(1,'Failed to define %s\n',variables{lvar});
		must_abort=1;
	    elseif strncmp(variables(lvar),'detsensefun',11)
		fprintf(1,'Failed to define %s\n',variables{lvar});
		must_abort=1;
	    end
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
%	fprintf(1,'Converted phasorsurface\n');
%	phasorsurface
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



fprintf('Allocating grid...');
fdtdgrid = initialisesplitgrid(I,J,K,Dxl,Dxu,Dyl,Dyu,Dzl,Dzu);
fprintf('Done\n');

%[epso , muo , c] = import_constants;

%We now implement an incident planewave using a total/scattered field
%formulation. For k<K1 we have scattered and for k>=K1 we have total.
%The field is coupled into system by means of update equations to_l maintain
%consistency.

%correct:
%refractive_index = sqrt(real(epsr(1)));

%omega_an = 2*pi*f_an;
%lambda_an = c/(f_an*refractive_index);
%wave_num_an = 2*pi/lambda_an;%wave number in m^-1

%fprintf('Initialising source field...');

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

%Set up the Ksource field. This has to be defined on a 2d array
%over the range (I0,I1)x(J0,J1). We calculate the field values
%assuming an origin for the illumination
i_source = (interface.I0(1):interface.I1(1)) - illorigin(1);
j_source = (interface.J0(1):interface.J1(1)) - illorigin(2);
k_source = interface.K0(1) - illorigin(3);



if interface.K0(2)

    %Ex, K0
    [x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ex');
    z = z + z_launch;
    ex_coords.x = x;
    ex_coords.y = y;
    ex_coords.z = z;

    %Ey, K0
    [x,y,z] = yeeposition(i_source,j_source,k_source,delta,'Ey');
    z = z + z_launch;
    ey_coords.x = x;
    ey_coords.y = y;
    ey_coords.z = z;
else
    ex_coords = [];
    ey_coords = [];
end

tvec_E = (1:Nt)*dt;

fvec_E = (0:(numel(tvec_E)-1))/numel(tvec_E)/diff(tvec_E(1:2));
fvec_E = fftshift(fvec_E);
ind0 = find(fvec_E==0);
fvec_E(1:(ind0-1)) = fvec_E(1:(ind0-1)) -fvec_E(ind0-1) - fvec_E(ind0+1);
fvec_E = ifftshift(fvec_E);

[epso , muo , c] = import_constants;
refractive_index = sqrt(real(epsr(1)));
lambda_an = c/(f_an*refractive_index);
hwhm = lambda_an^2/((c/refractive_index)*wavelengthwidth)*2*sqrt(log(2)/pi);
to_l = hwhm*sqrt(log(1e8)/pi);
