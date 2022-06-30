%function [verts,facets] = triangulateplane(I0,I1,J0,J1,skip,orientation)
function [verts,facets] = triangulateplane(I0,I1,J0,J1,skip,orientation)

%I0 = 1;I1 = 5;
%J0 = 2;J1 = 6;
%skip = 1;
%orientation = 1;
    if I1==I0 || J1==J0
	verts=[];facets=[];
    else
    ivec = I0:skip:I1;
    if ivec(end)~=I1
	ivec = [ivec I1];
    end
    
    jvec = J0:skip:J1;
    if jvec(end)~=J1 
	jvec = [jvec J1];
    end
    
    [VI,VJ] = ndgrid(ivec,jvec);
    
    verts = [VI(:) VJ(:)];

    facets = delaunay(verts(:,1),verts(:,2));
    %now check the orientation
    for ifac = 1:size(facets,1)
    %for ifac=1
	vec1 = verts(facets(ifac,2),:)-verts(facets(ifac,1),:);
	vec2 = verts(facets(ifac,3),:)-verts(facets(ifac,1),:);
	if (vec1(1)*vec2(2)-vec2(1)*vec1(2))*orientation < 0
	    ft = facets(ifac,:);
	    facets(ifac,2)=ft(3);
	    facets(ifac,3)=ft(2);
	end
	
    end

    
    verts = [verts(:,1) verts(:,2) zeros(size(verts(:,1)))];
    end
