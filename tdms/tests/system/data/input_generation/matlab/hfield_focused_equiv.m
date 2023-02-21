%function [E] = efield(X,Y,Z);
%
%This function is called by iteratefdtd_matrix to set the electric
%field source terms
%
function [E] = hfield_focused_equiv(X,Y,Z)


    [m,n] = size(X);

    for i=1:m
	for j=1:n
	    %[Etmp,d] = RichardsWolfZernike(pi/2,1,0,inf,0,0,1,10,0,0,X(i,j)*1e6,X(i,j)*1e6,Y(i,j)*1e6,Y(i,j)*1e6,Z*1e6,Z*1e6,1,1,1,10);
	    %the following -1* terms are there to account for the
            %fact that the fdtd code employs the exp(+jwt) convention
	    %E{1}(i,j) = exp(sqrt(-1)*k*Z);
	    %E{2}(i,j) = 0;
	    %[Etmp,Htmp,d] = NLayerIllumination(alpha,k,nvec,dbs,gbs,hvec,X(i,j),X(i,j),Y(i,j),Y(i,j),Z,Z,Nx,Ny,Nz,nint);

	    E{1}(i,j) = 0;%zeros(size(Htmp{1}));
	    E{2}(i,j) = 0;%zeros(size(Htmp{2}));
	    E{3}(i,j) = 0;%zeros(size(Htmp{3}));
	end
    end
