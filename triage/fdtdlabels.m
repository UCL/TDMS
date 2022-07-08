function [d] = fdtdlabels(fdtdin)
%function [d] = fdtdlabels(fdtdin)
%
%fdtdin - input structure to fdtd program
%
%d - axis labels for use in plotintensity etc
    
    d = {fdtdin.grid_labels.x_grid_labels,fdtdin.grid_labels.y_grid_labels,fdtdin.grid_labels.z_grid_labels};
 
