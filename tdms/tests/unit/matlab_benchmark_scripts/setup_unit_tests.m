if ~exist('matlab_data', 'dir')
    mkdir('matlab_data');
end

create_structure_array;
create_tdms_object_data;
create_bad_tdms_object_data;
