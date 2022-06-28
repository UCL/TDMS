%function [vertices,facets] = gen_mesh(x,y,z)
%
%Construct a mesh representing the surface of a cuboid with vertex
%points taken from the vector x, y and z. These vectors are assumed
%to be monotomically increasing or decreasing (or else the surface
%will not be a cuboid). It handles the degenerate case of a
%singleton vector by returning a triangulated plane.
%
%Only one singleton vector is permitted.
%
%Facets is based upon a 0 indexing scheme because this code is used
%predominantly in C. Thus for use in functions like trimesh, 1
%should be added to the facets.
%
%For example,
%
%Construct a mesh of a square, centred on the origin, ranging from
%x=-1 to x=1:
%
%x = linspace(-1,1,11);
%y = linspace(-1,1,11);
%z = 0;
%[vertices,facets] = gen_mesh(x,y,z);
%trimesh(double(facets)+1,vertices(:,1),vertices(:,2),vertices(:,3));
%
function [vertices,facets] = gen_mesh(x,y,z)
    
    I0 = 1;
    I1 = length(x);
    J0 = 1;
    J1 = length(y);
    K0 = 1;
    K1 = length(z);
    
    [vertices, facets] = mesh(I0,I1,J0,J1,K0,K1);
    vec_x = x(vertices(:,1)).';
    vec_y = y(vertices(:,2)).';
    vec_z = z(vertices(:,3)).';
    
    vec_x = reshape(vec_x,length(vec_x),1);
    vec_y = reshape(vec_y,length(vec_y),1);
    vec_z = reshape(vec_z,length(vec_z),1);
    
    
    vertices = [vec_x vec_y vec_z];

    
