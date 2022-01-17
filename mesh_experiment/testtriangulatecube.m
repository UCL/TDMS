I0 = 1;I1 = 10;
J0 = 2;J1 = 7;
K0 = 3;K1 = 8;

skip = 1;
orientation = 1;

[verts,facets] = triangulatecube(I0,I1,J0,J1,K0,K1,skip,orientation);
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
