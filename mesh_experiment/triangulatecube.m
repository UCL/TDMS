%Calculates a mesh of triangles with outward facing normals
%
%The normal is defined by taking the cross product between vectors
%formed with the 1st and 2nd vertices and the 1st and 3rd vertices.
%
%function [verts,facets] = triangulatecube(I0,I1,J0,J1,K0,K1,skip)
function [verts,facets] = triangulatecube(I0,I1,J0,J1,K0,K1,skip)

    Vyz = [0 0 -1;0 1 0;1 0 0];%I0->K0, I1->K1, 
    Fyz = [1 0 0;0 0 1;0 1 0];
    
    Vxz = [1 0 0;0 0 -1;0 1 0];%J0->K0, J1->K1, 
    Fxz = [1 0 0;0 0 1;0 1 0];
    
    [vertsxy1,facetsxy1] = triangulateplane(I0,I1,J0,J1,skip,-1);
    [vertsxy2,facetsxy2] = triangulateplane(I0,I1,J0,J1,skip, 1);
    if ~isempty(vertsxy1)
	vertsxy1(:,3) = K0;
	vertsxy2(:,3) = K1;
    end
    
    [vertsyz1,facetsyz1] = triangulateplane(K0,K1,J0,J1,skip,-1);
    [vertsyz2,facetsyz2] = triangulateplane(K0,K1,J0,J1,skip, 1);
    
    if ~isempty(vertsyz1)
	vertsyz1 = (Vyz*vertsyz1.').';
	facetsyz1 = (Fyz*facetsyz1.').';
	vertsyz2 = (Vyz*vertsyz2.').';
	facetsyz2 = (Fyz*facetsyz2.').';
        
	vertsyz1(:,1) = I0;
	vertsyz2(:,1) = I1;
    end
    
    [vertsxz1,facetsxz1] = triangulateplane(I0,I1,K0,K1,skip,-1);
    [vertsxz2,facetsxz2] = triangulateplane(I0,I1,K0,K1,skip, 1);
    
    if ~isempty(vertsxz1)
	vertsxz1 = (Vxz*vertsxz1.').';
	facetsxz1 = (Fxz*facetsxz1.').';
	vertsxz2 = (Vxz*vertsxz2.').';
	facetsxz2 = (Fxz*facetsxz2.').';
        
	vertsxz1(:,2) = J0;
	vertsxz2(:,2) = J1;
    end
    
    if I0==I1
	vertsyz2=[];facetsyz2=[];
    elseif J0==J1
	vertsxz2=[];facetsxz2=[];
    elseif K0==K1
	vertsxy2=[];facetsxy2=[];
    end
   
    verts = [vertsxy1;vertsxy2;vertsyz1;vertsyz2;vertsxz1;vertsxz2];
    facets = [facetsxy1;facetsxy2+size(vertsxy1,1);facetsyz1+size(vertsxy1,1)+size(vertsxy2,1);facetsyz2+size(vertsxy1,1)+size(vertsxy2,1)+size(vertsyz1,1);facetsxz1+size(vertsxy1,1)+size(vertsxy2,1)+size(vertsyz1,1)+size(vertsyz2,1);facetsxz2+size(vertsxy1,1)+size(vertsxy2,1)+size(vertsyz1,1)+size(vertsyz2,1)+size(vertsxz2,1)];
    
    %find any duplicate vertices
    %the starting point is that all duplicates occur on the edge,
    %thus they will have two extreme coordinates
    
    %inds= cell(1,12);% one vector for each edge
    edges = [I0 J0 0;I0 J1 0;I1 J0 0;I1 J1 0;I0 0 K0;I0 0 K1;I1 0 K0;I1 0 K1;0 J0 K0;0 J0 K1;0 J1 K0;0 J1 K1];
    dupes = [];
    for ie=1:size(edges,1)
	
	nz = find(edges(ie,:));
	%iz = find(~edges(ie,:));
	vl = find( (verts(:,nz(1))==edges(ie,nz(1))).*(verts(:,nz(2))==edges(ie,nz(2))) );
	dupes = [dupes; vl];
    end
    dupes=unique(dupes);
    %save dupes dupes;
    %now work out from this list which vertices are duplicated
    dupmap = zeros(numel(dupes),2);
    verts_c = verts;
    for ivl1 = 1:numel(dupes)
	for ivl2 = (ivl1+1):numel(dupes)
	    if isequal( verts(dupes(ivl1),:),verts_c(dupes(ivl2),:) )
		dupmap(ivl1,min(find(~dupmap(ivl1,:))))=dupes(ivl2);
		%if ivl1==1
		%    verts(dupes(ivl1),:)
		%    verts_c(dupes(ivl2),:)
		%    dupmap(ivl1,:)
		%end
		
		verts_c(dupes(ivl2),:)=[sqrt(-1) sqrt(-1) sqrt(-1)];
	    end
	end
    end
 
    %now work out the new mapping
    %in particular, each index  dupes(dupmap(i,:)) should be remapped to dupes(i)
    if 1
	mapping = zeros(numel(dupes),2);
	for id = 1:numel(dupes)
	    mapping(id,1)=dupes(id);
	    ind = find( dupmap == dupes(id) );
	    if isempty(ind)
		mapping(id,2)=mapping(id,1);
	    else
		[ii,jj] = ind2sub(size(dupmap),ind);
		mapping(id,2)=dupes(ii);
	    end
	end
    end
    %the second mapping takes account of the fact that we will
    %delete duplicated vertices
    
    vertsn = verts;
    for im=1:size(mapping,1)
	if mapping(im,1)~=mapping(im,2)
	    vertsn(mapping(im,1),:)=sqrt(-1);
	end
    end
    mapping2=zeros(numel(find(vertsn(:,1)~=sqrt(-1))),2);
    vertsf = [];
    for iv = 1:size(vertsn,1)
	if(vertsn(iv,1)~=sqrt(-1));
	    vertsf = [vertsf;vertsn(iv,:)];
	    mapping2(size(vertsf,1),1)=iv;
	    mapping2(size(vertsf,1),2)=size(vertsf,1); 
	end
	
    end

    %now do the remapping
    for im=1:size(mapping,1)
	if mapping(im,1)~=mapping(im,2)
	    %find all facets which index the duplicated vertex
	    facets(find(facets==mapping(im,1)))=mapping(im,2);
	end
    end
    
    for im=1:size(mapping2,1)
	if mapping2(im,1)~=mapping2(im,2)
	    %find all facets which index the duplicated vertex
	    facets(find(facets==mapping2(im,1)))=mapping2(im,2);
	end
    end
    verts = vertsf;
    
    %save dupmap dupmap dupes verts mapping vertsn vertsf mapping2;
