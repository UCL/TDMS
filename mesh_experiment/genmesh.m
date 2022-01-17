%function [vertices, facets] = genmesh([I0 I1 J0 J1 K0 K1])
function [vertices, facets] = genmesh(cuboidlims)
    
    if length(cuboidlims) ~= 6
	error('Failed to specify cuboid dimensions correctly');
    end
    
    [vertices, facets] = mesh(cuboidlims(1),cuboidlims(2),cuboidlims(3),cuboidlims(4),cuboidlims(5),cuboidlims(6));
	
