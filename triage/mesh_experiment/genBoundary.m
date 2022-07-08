function [verts,facets] = genBoundary(I0,I1,K0,K1)

%I0=2;
%I1=31;
%K0=3;
%K1=45;

Nverts = 2*(K1-K0) + 2*(I1-I0);
Nfacets = Nverts;

verts=zeros(Nverts,3);
facets=zeros(Nfacets,3);

verts_t=zeros(Nverts+4,3);
facets_t=zeros(Nfacets,3);
fcount=1;

tverts_i=(I0:I1);
tverts_k=(K0:K1);

nextv=1;Nv=I1-I0+1;
verts_t( nextv:(nextv+Nv-1),1) = I0:I1;
verts_t( nextv:(nextv+Nv-1),2) = 1;
verts_t( nextv:(nextv+Nv-1),3) = K0;

for ifind=1:(Nv-1)
	facets_t(fcount,1)=nextv+ifind-1;facets_t(fcount,2)=facets_t(fcount,1)+1;facets_t(fcount,3)=facets_t(fcount,1)+1;	
	fcount=fcount+1;
end

nextv=nextv+Nv;Nv=I1-I0+1;
verts_t( nextv:(nextv+Nv-1),1) = fliplr(I0:I1);
verts_t( nextv:(nextv+Nv-1),2) = 1;
verts_t( nextv:(nextv+Nv-1),3) = K1;
for ifind=1:(Nv-1)
	facets_t(fcount,1)=nextv+ifind-1;facets_t(fcount,2)=facets_t(fcount,1)+1;facets_t(fcount,3)=facets_t(fcount,1)+1;
	fcount=fcount+1;
end

nextv=nextv+Nv;Nv=K1-K0+1;
verts_t( nextv:(nextv+Nv-1),1) = I0;
verts_t( nextv:(nextv+Nv-1),2) = 1;
verts_t( nextv:(nextv+Nv-1),3) = fliplr(K0:K1);
for ifind=1:(Nv-1)
	facets_t(fcount,1)=nextv+ifind-1;facets_t(fcount,2)=facets_t(fcount,1)+1;facets_t(fcount,3)=facets_t(fcount,1)+1;
	fcount=fcount+1;
end

nextv=nextv+Nv;Nv=K1-K0+1;
verts_t( nextv:(nextv+Nv-1),1) = I1;
verts_t( nextv:(nextv+Nv-1),2) = 1;
verts_t( nextv:(nextv+Nv-1),3) = K0:K1;
for ifind=1:(Nv-1)
	facets_t(fcount,1)=nextv+ifind-1;facets_t(fcount,2)=facets_t(fcount,1)+1;facets_t(fcount,3)=facets_t(fcount,1)+1;
	fcount=fcount+1;
end

ind_map=ones(I1-I0+1,K1-K0+1)*-1;

vertc=1;
for fcount=1:Nfacets
    for vcount=1:2
	vind_i=verts_t(facets_t(fcount,vcount),1);
	vind_k=verts_t(facets_t(fcount,vcount),3);
	if ind_map(vind_i-I0+1,vind_k-K0+1)==-1
	    ind_map(vind_i-I0+1,vind_k-K0+1)=vertc;
	    verts(vertc,:)=[vind_i 1 vind_k];
	    vertc=vertc+1;
	end
	verti=ind_map(vind_i-I0+1,vind_k-K0+1);
	facets(fcount,vcount)=verti;
    end
    facets(fcount,3)=verti;
end


