%function [material_matrix,composition_matrix] = read_material_data(filename)
%
%Can open an ascii file with the following specifications or a mat
%file which is described after. The functions determines which is
%which based upon the file extension. A file extension of .mat is
%assumed to be a mat file.
%
%================================================================================
%Ascii file:
%================================================================================
%
%Opens 'filename' and reads the grid material data
%
%material_matrix is a matrix with rows of the form:
%           [material_identifier eps_r mu_r sigma_x .... sigma_z*]
%
%composition_matrix is a matrix of the form:
%
%           [i j k material_identifier]
%
%The expected file format is as follows:
%
%-all lines commencing with % are ignored
%
%-the file is separated into regions each being defined by a heading. The headings
% currently in use are:
%                       materials
%                       cells
%
%-a material is defined using the following format:
%   <material identifier> <eps_r> <mu_r> <sigma_x> <sigma_y> <sigma_z> <sigma_x*> <sigma_y*> <sigma_z*> 
%
%   or
%
%   <material identifier> <eps_r> <mu_r> <vc> <wp> <sigma_x> <sigma_y> <sigma_z> <sigma_x*> <sigma_y*> <sigma_z*> 
%
%   where a dispersive material is to be defined.
%
%   For example, the following entry represents free space:
%   0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
%
%   Free space is by default material 0 and need not be entered. It may be entered however.
%
%   These entries should appear after the materials heading and before any other heading.
%
%-the contents of a cell is specified by a line with the following format:
%   i j k <material identifier>
%
%   the indices i j k are relative to an origin which is defined as the first non PML
%   cell in the grid. The first means that with each index the minimum such that the 
%   cell is non-PML. This means that the origin has global index (Dxl, Dyl, Dzl).
%
%   indexing should start from (1,1,1) s.t. the cell (1,1,1) is the first non-PML
%   cell.
%
%   These entries should appear after the cells heading and before
%   any other heading.
%
%================================================================================
%Mat file
%================================================================================
%
%The mat file should contain 2 matrices, material_matrix and
%composition_matrix which are essentially matrix versions of what
%is oulined above.
function [material_matrix,composition_matrix] = read_material_data(filename)

%first determine the type of file we have to open
file_mat   =  0;
file_ascii =  1;
file_type  = -1;
if length(filename) > 4
    if strncmp(filename((length(filename)-3):length(filename)),'.mat',4)
	file_type = file_mat;
    else
	file_type = file_ascii;
    end
else
    file_type = file_ascii;
end

if file_type == file_ascii
    
			       
    %but need to store whether or not any
    %entries are actually made
    t_composition_filled = 0;
    t_material_filled = 0;
    
    
    fid = fopen(filename);
    
    if fid == -1
	error(sprintf('Failed to open %s\n',filename));
    end
    
    buffer = fgets(fid);
    
    material_counter = 1;
    composition_counter = 1;
    
    t_material = zeros(1,11);
    t_composition = zeros(1,4);%used to build the above matrices - used
			       %for efficiency
        			       
    mode = 0; %1 for materials, 2 for cells
    
    while buffer ~= -1 %while not at end of file
	if ~strncmp(buffer,'%',1) & ~is_white_space(buffer)%if not a comment
	    if strncmp(buffer,'materials',length('materials'))
		mode = 1;
		buffer = fgets(fid);
	    elseif strncmp(buffer,'cells',length('cells'))
		mode = 2;
		buffer = fgets(fid);
	    end
	    
	    if buffer ~= -1 & ~strncmp(buffer,'%',1) & ~is_white_space(buffer)
		if mode==1
		    [material_data,count] = sscanf(buffer,'%d %e %e %e %e %e %e %e %e %e %e');
		    if count ~= 9 & count ~= 11
			error(sprintf('Encountered incorrect line in %s - %s',filename,buffer));    
		    end
		    if count==11
			t_material(material_counter,:) = material_data.';
		    else
			t_material(material_counter,1:3)  = material_data(1:3).';
			t_material(material_counter,4:5)  = [0 0];
			t_material(material_counter,6:11) = material_data(4:9).';
		    end
		    t_material_filled = 1;
		    material_counter = material_counter + 1;   
		elseif mode==2
		    [composition_data,count] = sscanf(buffer,'%d %d %d %d');
		    if count ~= 4
			error(sprintf('Encountered incorrect line in %s - %s',filename,buffer));        
		    end
		    t_composition(composition_counter,:) = composition_data.';
		    t_composition_filled = 1;
		    composition_counter = composition_counter + 1;
		end
	    end
	end
	if buffer ~= -1
	    buffer = fgets(fid);      
	end
    end
else
    dat = load(filename);
    t_material = dat.material_matrix;
    %now update material_matrix to ensure it has 11 columns
    [nrows, ncols] = size(t_material);
    if ncols==9
        t_material = [t_material(:,1:3) zeros(nrows,2) t_material(:,4:9)];
    end     
    t_composition = dat.composition_matrix; 
        
    t_material_filled = 1;
    t_composition_filled = 1;
end


%now just check that every material identifier used in the cells section has been defined in the materials section
[m,n] = size(t_composition);
for i=1:m
    if sum( t_material(:,1) == t_composition(i,4) ) == 0
        error(sprintf('Material %d does not exist',t_composition(i,4)));
    end
end

if t_material_filled
    material_matrix = zeros(size(t_material));
    material_matrix(:,:,:,:,:,:,:,:,:,:,:) = t_material(:,:,:,:,:,:,:,:,:,:,:);
    clear t_material;
else
    material_matrix =[];
end

if t_composition_filled
    composition_matrix = zeros(size(t_composition));
    composition_matrix(:,:,:,:) = t_composition(:,:,:,:);
    clear t_composition;
else
    composition_matrix = [];
end

