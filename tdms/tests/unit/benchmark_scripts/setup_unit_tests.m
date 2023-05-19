if ~exist('unit_test_data', 'dir')
    mkdir('unit_test_data');
end

create_structure_array;
create_class_data;
create_bad_class_data;

exit;
