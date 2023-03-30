function create_tdfield_data(output_name, I_tot, J_tot, Nt)
    %% Creates valid time-domain field data and saves it to a .mat folder.
    % The field created is just an array of 0s of the correct size [I_tot+1, J_tot+1, Nt].

    % The input file we are using will accept td-fields of dimensions
    % 277-1-500 as valid inputs, which are the defaults set here.
    if ~exist('I_tot', 'var')
        I_tot = 276;
    end
    if ~exist('J_tot', 'var')
        J_tot = 0;
    end
    if ~exist('Nt', 'var')
        Nt = 500;
    end

    exi = zeros(I_tot + 1, J_tot + 1, Nt);
    eyi = zeros(I_tot + 1, J_tot + 1, Nt);
    save(output_name, 'exi', 'eyi');
end
