I0 = 1;I1 = 5;
J0 = 2;J1 = 6;
skip = 2;
orientation = 1;

[verts,facets] = triangulateplane(I0,I1,J0,J1,skip,orientation);

if 0
    Vyz = [0 0 -1;0 1 0;1 0 0];%I0->K0, I1->K1, 
    Fyz = [1 0 0;0 0 1;0 1 0];
    verts = (Vyz*verts.').';
    facets= (Fyz*facets.').';
end
Vxz = [1 0 0;0 0 -1;0 1 0];%J0->K0, J1->K1, 
Fxz = [1 0 0;0 0 1;0 1 0];
verts = (Vxz*verts.').';
facets= (Fxz*facets.').';

figure(1);clf;
trisurf(facets,verts(:,1),verts(:,2),verts(:,3),zeros(size(verts(:,1))));

for ifac = 1:size(facets,1)
    vec1 = verts(facets(ifac,2),:)-verts(facets(ifac,1),:);
    vec2 = verts(facets(ifac,3),:)-verts(facets(ifac,1),:);
%    k = vec1(1)*vec2(2)-vec2(1)*vec1(2);
kvec = cross(vec1,vec2);
    centroid = mean(verts(facets(ifac,:),:),1);
    
    line([centroid(1) centroid(1)+kvec(1)],[centroid(2) centroid(2)+kvec(2)],[centroid(3) centroid(3)+kvec(3)]);

    
end
hx = xlabel('x');
hy = ylabel('y');
hz = zlabel('z');
