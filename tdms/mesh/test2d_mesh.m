I0=1;
I1=5;
K0=2;
K1=6;

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

theta_sum=0;centI=(I1+I0)/2;centK=(K0+K1)/2;
vert3=[centI 1 centK];
figure(1);clf;
for faci=1:size(facets,1)
    vert1=verts(facets(faci,1),:);
    vert2=verts(facets(faci,2),:);
    
    vr1=vert1-vert3;
    vr2=vert2-vert3;
    
    r1=sqrt(sum( (vert1-vert3).^2));;
    r2=sqrt(sum( (vert2-vert3).^2));;
    e=sqrt(sum( (vert2-vert1).^2));;
    
    theta=acos(-(e^2-r1^2-r2^2)/2/r1/r2);
    theta_sum=theta+theta_sum;
    
    hl=line([vert1(1) vert2(1)],[vert1(2) vert2(2)],[vert1(3) vert2(3)]);hold on;
    plot3(vert1(1),vert1(2),vert1(3),'.');
    plot3(vert2(1),vert2(2),vert2(3),'.');
    
    
end
 view([0 -1 0]);


%verts( ((I1-I0)+1),1) = (I0:I1)-1;

