function create_gridfile(output_name)
    %% Create a valid gridfile to pass into iteratefdtd_matrix.
    % The gridfile will be saved under the output_name provided.

    composition_matrix = [];
    material_matrix = [1 0 1 0 0 0 0 0 0 0 0];
    save(output_name, 'composition_matrix', 'material_matrix');
end
