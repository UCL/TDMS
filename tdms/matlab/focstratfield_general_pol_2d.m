%function [Ep,Em] = focstratfield_general_pol(vertices,nvec,hvec,NA,lambda,ntheta,nphi,polfun)
%
%Inputs:
%vertices is an N x 3 matrix of positions at which the field is to
%         be calculated
%
%nvec     is a vector of the form [n_0,n_1,...,n_N] defining the
%         refractive indices of the layers of the stratified medium
%
%hvec     is a vector of the form [h_1,h_2,...,h_N] defining the
%         z locations of the unrotated stratified medium
%
%NA       The numerical aperture of the lens in medium n_0. If NA
%         is a vector it must have two elements. In this case the
%         first element specifies the inside of a transmitting
%         annulus and the second element specifies the outer extent
%         of the annulus.
%
%lambda   wavelength
%
%ntheta   number of angles in theta to integrate over
%
%nphi     number of angles in phi to integrate over
%         if nphi is set to [], it is assumed that a cylindrical
%         lens is consisdered and so only phi=0 and phi=pi is considered.
%
%Output:
%Ep       N x 3 matrix of forward propagating fields at each vertex
%
%Em       N x 3 matrix of backward propagating fields at each
%         vertex
%
%Based on the theory published in ...
function [Ep,Em] = focstratfield_general_pol(vertices,nvec,hvec,NA,lambda,ntheta,nphi,polfun)

%To do
%1) Update reference to published paper above
%2) Update equation numbers cited in comments below
    
%The following parameters are currently not set as the features
%they provide have not been fully tested

%theta_sm rotation angle of the stratified medium about z=zr
theta_sm=0;

%zr       point on the optical axis about which stratified medium
%         is rotated
zr=0;
%end of untested parameters

%perform basic error checking
%nvec must have one more element than hvec
if numel(nvec) ~= numel(hvec)+1
    error('nvec must have one more element than hvec');
end

%now check that there are no repeated entries in hvec, remove
%them if there are
[hvec_dash,ia,ic] = unique(hvec);

if numel(hvec_dash) ~= numel(hvec)
    %there were some duplicate entries in hvec
    error('hvec should contain unique values');
end

%hvec must be in ascending order
if ~all(hvec==sort(hvec))
    error('hvec must be in ascending order');
end

%check if NA is valid
if numel(NA)~=1
    if numel(NA)~=2
	error('NA must have either 1 or 2 entries');
    else
	if ~all(NA/nvec(1)<1)
	    error('NA/nvec(1) must be < 1');
	end
	
	if NA(1)>=NA(2)
	    error('NA must be in ascending order');
	end
    end
else
    if NA<0
	error('NA must be > 0');
    end
    
    if NA/nvec(1)>1
	error('NA/nvec(1) must be < 1');
    end
end
%end of error checking

%Determine the number of unique values of z
[zplanelocs,iunique,imaptozloc] = unique(vertices(:,3));
Nplanes = numel(zplanelocs);

%now we generate the theta and phi vectors, in particular:
%
%phi_vec: vector of azimuthal angles considered in angular
%spectrum of focussed beam
%
%theta_vec: vector of polar angles, measured between the ray
%and z-axis, of plane waves considered in the angular
%spectrum of focussed beam
if isempty(nphi)
    %[phi_vec1,phi_w1]=gauss_legendre(-1e-9,1e-9,10);
    %[phi_vec2,phi_w2]=gauss_legendre(pi+-1e-9,pi+1e-9,10);
    %phi_vec = [phi_vec1 phi_vec2];
    %phi_w = [phi_w1 phi_w2];
    phi_vec = [0 pi];
    phi_w = [1 1]*1/2;
else
    [phi_vec,phi_w]=gauss_legendre(0,2*pi,nphi);
end

if numel(NA)==2
    [theta_vec,theta_w]=gauss_legendre(asin(NA(1)/nvec(1)),asin(NA(2)/nvec(1)),ntheta);
else
    [theta_vec,theta_w]=gauss_legendre(0,asin(NA/nvec(1)),ntheta);
end

%Create matrices allowing for integration over theta and phi later
[PHI,THE] = ndgrid(phi_vec,theta_vec);
[WPHI,WTHE] = ndgrid(phi_w,theta_w);

%convert all to vectors as we will need to combine them later
%with the observation coordinates
%
%This is for computational convenience
WPHI = WPHI(:);
WTHE = WTHE(:);
PHI = PHI(:);
THE = THE(:);

%generate the propagation vectors (N x 3)
%
%S is thus a matrix of dimensions N x 3 where each column is a
%unit vector describing the direction of propagation of a plane
%wave in the focal region. See the definition of
%$\boldsymbol{s}$ just after Eq (3).
S = [cos(PHI).*sin(THE) sin(PHI).*sin(THE) cos(THE)];

%generate polarisation vectors - this needs to be modified to
%model other beam polarisations
%E = (sqrt(cos(THE)).*sin(THE)*ones(1,3)).*[ 0.5*((1+cos(THE))+cos(2*PHI).*(cos(THE)-1)) 0.5*(sin(2*PHI).*(cos(THE)-1)) -sin(THE).*cos(PHI)];

%This applies specifically to an x-polarised plane wave being
%focussed. This is Eq (5) but without the $\sqrt{\cos\theta}$
%term
Ex = [ 0.5*((1+cos(THE))+cos(2*PHI).*(cos(THE)-1)) 0.5*(sin(2*PHI).*(cos(THE)-1)) -sin(THE).*cos(PHI)];
Ey = [ 0.5*(sin(2*PHI).*(cos(THE)-1)) 0.5*((1+cos(THE))-cos(2*PHI).*(cos(THE)-1)) -sin(THE).*sin(PHI)];

%Calculate weights
[wx,wy] = polfun(THE,PHI);

E = Ex.*(wx*[1 1 1]) + Ey.*(wy*[1 1 1]);

%Now we apply the  $\sqrt{\cos\theta}$ term from Eq (5), the
%$\sin\theta$ term from Eq (4) and the differentials $d\theta$
%and $d\phi$.
if isempty(nphi)
    WGAUSS = (sqrt(cos(THE)));%.* WGAUSS;
else
    WGAUSS = (sqrt(cos(THE)).*sin(THE));%.* WGAUSS;
end
%transformation matrix used to model rotation of medium, Eq(1)
%except the rotation is being applied in the opposite
%direction, hence opposing signs in the sine terms
T = [cos(theta_sm) 0 sin(theta_sm)
     0           1 0
     -sin(theta_sm) 0 cos(theta_sm)];

%inverse of T
Ti = T';

%We now transform S and observation points into rotated coordinate
%system
St = (Ti*S')';
Et = (Ti*E')';
Rt = (Ti*vertices')';

%calculate the location of the  original origin in the
%reference frame O4 (see text for description of the different
%refrence frames), Eq (2) 
origin = (Ti*([0 0 0]-zr*1*[0 0 1])')' + zr*1*[0 0 1];

%The layers have the same z positions in reference frame O4, we
%want to know what they are in reference frame O5 which means
%subtracting the z-component of the old origin. See text just
%after Eq (2).
hvec = hvec - origin(3);

%Original vertices in order of increasing z
[z_ordered,ind_z_ordered] = sort(vertices(:,3));
R = vertices(ind_z_ordered,:);

%we want to make the algorithm blind to whether or not there is
%a stratified medium. So the first thing to do is create
%matching hvec and nvec that contain the first and final
%vertices.

%if the first vertex appears before the first interface, create
%a phantom interface
if isempty(hvec)
    if (R(end,3)-R(1,3))<eps
	hvec = R(end,3);
	nvec = [nvec nvec];
    else
	hvec = [R(1,3) R(end,3)];
	nvec = [nvec nvec nvec];
    end
else
    if (R(1,3) - hvec(1)) <= -eps
	hvec = [R(1,3) hvec];
	nvec = [nvec(1) nvec];
    end
    %We no longer do this in case there is high attenuation after the
    %final interface
    %	if (R(end,3) - hvec(end)) > eps;
    %	    hvec = [hvec R(end,3)];
    %	    nvec = [nvec nvec(end)];
    %	end
end

%Incident propagation vector at first interface
St_0_p = [St(:,1) St(:,2)  St(:,3)];    

%Reflected propagation vector at first interface
St_0_m = [St(:,1) St(:,2) -St(:,3)];    

%Transmitted  propagation vector after final interface
St_N_p = [nvec(1)/nvec(end)*St(:,1) nvec(1)/nvec(end)*St(:,2) sqrt( 1-(nvec(1)/nvec(end)*St(:,1)).^2 -(nvec(1)/nvec(end)*St(:,2)).^2)];

%Now we work out the unit vectors for the TE and TM electric field components.
%subscript 1 means incident
%Eq (6), additional code ensures the vectors have unit
%magnitude

%Field vectors before the first interface
TEhat_0_p = [-St(:,2) St(:,1) zeros(size(St(:,2)))];
inds = find(sqrt(sum(abs(TEhat_0_p.').^2)).'<eps);
TEhat_0_p(inds,1)=0;TEhat_0_p(inds,2)=1;TEhat_0_p(inds,3)=0;
TEhat_0_p_norm = sqrt(sum(abs(TEhat_0_p.').^2)).';
TEhat_0_p = TEhat_0_p./(TEhat_0_p_norm*ones(1,3));
%Calculate the component of the field parallel to the TE unit vector
TEmag_0_p = sum((TEhat_0_p.*E).').';

%Eq (7)
TMhat_0_p = cross(TEhat_0_p,St);
%Calculate the component of the field parallel to the TM unit vector
TMmag_0_p = sum((TMhat_0_p.*E).').';

%calculate the angle of incidence
THEi = acos(St(:,3));

%The general strategy for the calculation of RTE, RTM, TTE, TTM
%and the interpolated versions is to calculate the transmission
%coefficient just after the final physical interface. These
%coefficients are called:
%
%TTEi_final_boundary
%TTMi_final_boundary
%
%We then work through each value of z where the field is to be
%evaluated and we form the matrix S which links the positive
%and negative going fields at this plane, with the field
%transmitted just after the final interface. So the positive
%and negative going fields are always evaluated on the negative-z
%side of the interface.
RTE = zeros(numel(iunique),numel(THEi));
%Same as above but for TM
RTM = zeros(numel(iunique),numel(THEi));
%Same as above but transmission coefficient for TE
TTE = ones(numel(iunique),numel(THEi));
%Same as above but for TM
TTM = ones(numel(iunique),numel(THEi));

%refractive index just before the ith interface. This will be
%used to calculate the propagation vector or planes waves
%incident upon that interface
nvec_internal = zeros(1,numel(iunique));

Ep = zeros(size(vertices));Em = zeros(size(vertices));

%find the unique value of zplanelocs which is either before or
%at the final interface
ifinalz = max(find(zplanelocs<=hvec(end)));
if isempty(ifinalz)
    ifinalz=0;
end

if ~isempty(hvec)
    %The point of this is that we don't want to evaluate the
    %transmission and reflection coefficients at all angles in
    %THEi, since many are repeated.
    %
    %Instead, we calculate the reflection and transmission
    %coefficients across a range of theta angles covering the
    %range of THEi, then we interpolate for each value in THEi afterwards.
    %
    %If we didn't allow for rotating the stratified medium we
    %wouldn't need to perform this interpolation.
    %theta_i = linspace(min(THEi(:)),max(THEi(:)),ntheta);
    %Do this if interpolation isn't necessary
    theta_i = theta_vec;
    RTEi = zeros(numel(iunique),numel(theta_i));
    TTEi = zeros(numel(iunique),numel(theta_i));
    RTMi = zeros(numel(iunique),numel(theta_i));
    TTMi = zeros(numel(iunique),numel(theta_i));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %transmission just after the final interface
    TTEi_final_boundary = zeros(numel(theta_i));
    TTMi_final_boundary = zeros(numel(theta_i));
    for it=1:numel(theta_i)
	%These are the reflectance and transmission
	%coefficients at the first and last interfaces
	[STE] = multi_layer( hvec, nvec, lambda, theta_i(it), 'TE');
	TTEi_final_boundary(it) = 1/STE(1,1);
		
	[STM] = multi_layer( hvec, nvec, lambda, theta_i(it), 'TM');
	TTMi_final_boundary(it) = 1/STM(1,1);
		
	%nvec has been adjusted by this stage in the code to include a
	%phantom interface, if necessary, to accomodate an observation
	%point before the first interface. The final element of
	%nvec_internal should be the refractive index before
	%the first interface, either phantom or physical.
	
	%now we need to iterate through all possible
	%z-locations
	for iz=1:ifinalz
	    %z location of current vertices
	    zloc = vertices(iunique(iz),3);
	    %first find if zloc is equal to any elements within
	    %hvec
	    ind_zloc = find(abs(zloc - hvec)<eps);
	    if ~isempty(ind_zloc)
		hvec_temp = hvec(ind_zloc:end);
		nvec_temp = nvec(ind_zloc:end);
	    else
		ind_zloc = min(find(zloc < hvec));
		hvec_temp = [zloc hvec(ind_zloc:end)];
		nvec_temp = [nvec(ind_zloc) nvec(ind_zloc:end)];
	    end
	    if it==1
		nvec_internal(iz) = nvec_temp(1);
	    end
	    
	    theta_layer = asin(nvec(1)*sin(theta_i(it))/nvec_temp(1));
	    
	    [STE] = multi_layer( hvec_temp, nvec_temp, lambda, theta_layer, 'TE');
	    RTEi(iz,it) = STE(2,1)*TTEi_final_boundary(it);
	    TTEi(iz,it) = STE(1,1)*TTEi_final_boundary(it);
	    
	    [STM] = multi_layer( hvec_temp, nvec_temp, lambda, theta_layer, 'TM');
	    RTMi(iz,it) = STM(2,1)*TTMi_final_boundary(it);
	    TTMi(iz,it) = STM(1,1)*TTMi_final_boundary(it);
	    
	end
	for iz=(ifinalz+1):numel(iunique)
	    %z location of current vertices
	    zloc = vertices(iunique(iz),3);
	    
	    hvec_temp = [hvec(end) zloc];
	    nvec_temp = [nvec(end) nvec(end) nvec(end)];
	    if it==1
		nvec_internal(iz) = nvec_temp(1);
	    end
	    
	    theta_layer = asin(nvec(1)*sin(theta_i(it))/nvec_temp(1));
	    
	    [STE] = multi_layer( hvec_temp, nvec_temp, lambda, theta_layer, 'TE');
	    RTEi(iz,it) = 0;
	    TTEi(iz,it) = TTEi_final_boundary(it)/STE(1,1);
	    
	    [STM] = multi_layer( hvec_temp, nvec_temp, lambda, theta_layer, 'TM');
	    RTMi(iz,it) = 0;
	    TTMi(iz,it) = TTMi_final_boundary(it)/STM(1,1);
	    
	end
    end
    %Perform the inerpolation. Currently, each element of THEi
    %appears in theta_i so no values are actually interpolated.
    for iz=1:numel(iunique)
	RTE(iz,:) = interp1(theta_i,squeeze(RTEi(iz,:)),THEi(:),'pchip');
	TTE(iz,:) = interp1(theta_i,squeeze(TTEi(iz,:)),THEi(:),'pchip');
	RTM(iz,:) = interp1(theta_i,squeeze(RTMi(iz,:)),THEi(:),'pchip');
	TTM(iz,:) = interp1(theta_i,squeeze(TTMi(iz,:)),THEi(:),'pchip');
    end
end

%in some cases the (x,y) coordinates of the verties are
%repeated for different values of z. Here we check whether this
%is the case, and if so, we can save come computation below.

%first check whether there are the same number of vertices for
%each value of iz
Nix = zeros(1,numel(iunique));
for iz=1:numel(iunique)
    Nix(iz)=numel(find(imaptozloc.'==iz));
end
%Now check whether each plane of vertices contains the same
%transverse coordinates
transalign = 0;
if all(Nix==Nix(1))
    %now we check that all (x,y) coordinates are the same
    transalign = 1;
    for iz=1:numel(iunique)
	if transalign
	    if iz==1
		ixvec_0=find(imaptozloc.'==iz);
		if ~all(ixvec_0==(1:numel(ixvec_0)))
		    transalign = 0;
		end
	    else
		ixvec_1=find(imaptozloc.'==iz);
		if sqrt(sum(sum(abs(Rt(ixvec_1,1:2) - Rt(ixvec_0,1:2)).^2))) > eps
		    transalign = 0;
		end
	    end
	end
    end
end
%We have all information regarding planes waves in the spectrum
%and their reflection and transmission coefficients. Now we
%evaluate  Eq (3) for each vertex of interest
if transalign
    ixvec_all = zeros(numel(iunique),numel(ixvec_0));
    for iz=1:numel(iunique)
	ixvec_all(iz,:) = find(imaptozloc.'==iz);
    end
    for ix=1:numel(ixvec_0)
	PHASE1 = exp(sqrt(-1)*nvec(1)*2*pi/lambda*(Rt(ix,1:2)*St_0_p(:,1:2).' + hvec(1)*St_0_p(:,3).')).';
	for iz=1:numel(iunique)
	    %ixvec_1=find(imaptozloc.'==iz);
	    ixvec_1 = ixvec_all(iz,:);
	    computevecs = 1;
	    if iz>1
		if abs( nvec_internal(iz) - nvec_internal(iz-1) ) < eps
		    computevecs = 0;
		end
	    end
	    if computevecs
		TEhat_2 = TEhat_0_p;
		St_iz_p = [nvec(1)/nvec_internal(iz)*St(:,1) nvec(1)/nvec_internal(iz)*St(:,2) sqrt( 1-(nvec(1)/nvec_internal(iz)*St(:,1)).^2 -(nvec(1)/nvec_internal(iz)*St(:,2)).^2)];
		TMhat_2_p = cross(TEhat_0_p,St_iz_p);
		
		St_iz_m = [nvec(1)/nvec_internal(iz)*St(:,1) nvec(1)/nvec_internal(iz)*St(:,2) -sqrt( 1-(nvec(1)/nvec_internal(iz)*St(:,1)).^2 -(nvec(1)/nvec_internal(iz)*St(:,2)).^2)];
		TMhat_2_m = cross(TEhat_0_p,St_iz_m);
	    end
	    tmult_1 = WGAUSS.*WPHI.*WTHE.*PHASE1;
	    %perform the integration over the TE and TM
	    %components. 
	    
	    Ep(ixvec_1(ix),:) = (tmult_1.*TEmag_0_p.*TTE(iz,:).').'*TEhat_2 + (tmult_1.*TMmag_0_p.*TTM(iz,:).').'*TMhat_2_p;
	    Em(ixvec_1(ix),:) = (tmult_1.*TEmag_0_p.*RTE(iz,:).').'*TEhat_2 + (tmult_1.*TMmag_0_p.*RTM(iz,:).').'*TMhat_2_m;
	end
    end
else
    for iz=1:numel(iunique)
	for ix=find(imaptozloc.'==iz)
	    PHASE1 = exp(sqrt(-1)*nvec(1)*2*pi/lambda*(Rt(ix,1:2)*St_0_p(:,1:2).' + hvec(1)*St_0_p(:,3).')).';
	    tmult_1 = WGAUSS.*WPHI.*WTHE.*PHASE1;
	    %perform the integration over the TE and TM
	    %components. 
	    TEhat_2 = TEhat_0_p;
	    St_iz_p = [nvec(1)/nvec_internal(iz)*St(:,1) nvec(1)/nvec_internal(iz)*St(:,2) sqrt( 1-(nvec(1)/nvec_internal(iz)*St(:,1)).^2 -(nvec(1)/nvec_internal(iz)*St(:,2)).^2)];
	    TMhat_2 = cross(TEhat_0_p,St_iz_p);
	    Ep(ix,:) = (tmult_1.*TEmag_0_p.*TTE(iz,:).').'*TEhat_2 + (tmult_1.*TMmag_0_p.*TTM(iz,:).').'*TMhat_2;

	    St_iz_p = [nvec(1)/nvec_internal(iz)*St(:,1) nvec(1)/nvec_internal(iz)*St(:,2) -sqrt( 1-(nvec(1)/nvec_internal(iz)*St(:,1)).^2 -(nvec(1)/nvec_internal(iz)*St(:,2)).^2)];
	    TMhat_2 = cross(TEhat_0_p,St_iz_p);
	    Em(ix,:) = (tmult_1.*TEmag_0_p.*RTE(iz,:).').'*TEhat_2 + (tmult_1.*TMmag_0_p.*RTM(iz,:).').'*TMhat_2;
	end
    end
end

%convert back to the original coordinate system
Ep = -sqrt(-1)*nvec(1)/lambda*(T*Ep.').';Em = -sqrt(-1)*nvec(1)/lambda*(T*Em.').';

save fpgvars;
