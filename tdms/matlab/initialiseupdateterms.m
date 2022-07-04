%function [sigma, C, D, freespace] = initialiseupdateterms(R0, I,
%J, K, Dxl, Dxu, Dyl, Dyu, Dzl, Dzu, n, delta, dt, epsr_star, mur,
%multilayer, omega)
%
%Initialises the update equations for the Yee algorithm. Primary
%purpose is to implement the PML.
%
%Inputs - 
%
%R0 - maximum reflection tolerated
%I - the number of Yee cells (excluding PML) in the x direction
%J - the number of Yee cells (excluding PML) in the y direction
%K - the number of Yee cells (excluding PML) in the z direction
%Dx - Number of cells in x direction in a single PML layer
%Dy - Number of cells in y direction in a single PML layer
%Dz - Number of cells in z direction in a single PML layer
%n - order of increase in PML conductivity
%delta - spatial difference between grid points
%
function [sigma, C, D, freespace, conductive_aux, dispersive_aux] = initialiseupdateterms(R0, I, J, K, Dxl, Dxu, Dyl, Dyu, Dzl, Dzu, n, delta, dt, epsr_star, mur, multilayer, omega, kappa_max, vc_vec_i, wp_vec_i)
%    fprintf(1,'Entering initialiseupdateterms');
%define some constants
[epso , muo , c] = import_constants;    

%Initialise the sigma matrices
sigma_x.E = zeros(I+Dxl+Dxu+1,K+Dzl+Dzu+1);
sigma_x.H = zeros(I+Dxl+Dxu+1,K+Dzl+Dzu+1);
sigma_y.E = zeros(J+Dyl+Dyu+1,K+Dzl+Dzu+1);
sigma_y.H = zeros(J+Dyl+Dyu+1,K+Dzl+Dzu+1);
sigma_z.E = zeros(1,K+Dzl+Dzu+1);
sigma_z.H = zeros(1,K+Dzl+Dzu+1);


%process epsr
epsr_vec = real(epsr_star);
sigma_vec = (epsr_vec - epsr_star)*(sqrt(-1)*omega*epso);

epsr_z = zeros(1,K+Dzl+Dzu+1);
epsr_x = zeros(I+Dxl+Dxu+1,K+Dzl+Dzu+1);
epsr_y = zeros(J+Dyl+Dyu+1,K+Dzl+Dzu+1);

vc_vec = zeros(1,K+Dzl+Dzu+1);
wp_vec = zeros(1,K+Dzl+Dzu+1);

kappa_max_z = ones(1,K+Dzl+Dzu+1);

sigma_int_x = ones(size(sigma_x.E));
sigma_int_y = ones(size(sigma_y.E));
sigma_int_z = ones(size(sigma_z.E));



if isempty(multilayer)
    epsr = epsr_vec;
    epsr_z(:) = epsr;
    epsr_x(:) = epsr;
    epsr_y(:) = epsr;
    kappa_max_z(:) = kappa_max;
    
    vc_vec(:) = vc_vec_i;
    wp_vec(:) = wp_vec_i;
    
    sigma_int = sigma_vec;
    sigma_int_x(:) = sigma_int;
    sigma_int_y(:) = sigma_int;
    sigma_int_z(:) = sigma_int;
else    
    %epsr = 1;%this needs to be addressed
    epsr = epsr_vec(1);%this is less incorrect...
    start_ind = 1+Dzl;
    
    padded_multilayer = [multilayer (K+1)]+Dzl;
    
    for tc=1:length(padded_multilayer)
	for td=start_ind:(padded_multilayer(tc)-1)
%	    fprintf(1,'td: %d ',td);
	    epsr_z(td) = epsr_vec(tc);
	    vc_vec(td) = vc_vec_i(tc);
	    wp_vec(td) = wp_vec_i(tc);
	    kappa_max_z(td) = kappa_max(tc);
	    sigma_int_z(td) = sigma_vec(tc);
	end
%	fprintf(1,'\n');
	start_ind = padded_multilayer(tc) ;
    end

    for td=1:(Dzl)
%	fprintf(1,'td: %d ',td);
	epsr_z(td) = epsr_z(Dzl+1);
	kappa_max_z(td) = kappa_max_z(Dzl+1);
	sigma_int_z(td) = sigma_int_z(Dzl+1);
	vc_vec(td) = vc_vec(Dzl+1);
	wp_vec(td) = wp_vec(Dzl+1);
    end
%    fprintf(1,'\n');
    for td=(Dzl+K+1):(Dzl+K+1+Dzu)
	sigma_int_z(td) = sigma_int_z(Dzl+K);
    	epsr_z(td) = epsr_z(Dzl+K);
	kappa_max_z(td) = kappa_max_z(Dzl+K);
	vc_vec(td) = vc_vec(Dzl+K);
	wp_vec(td) = wp_vec(Dzl+K);
%	fprintf(1,'td: %d ',td);
    end
    sigma_int_x = ones(I+Dxl+Dxu+1,1)*sigma_int_z;
    sigma_int_y = ones(J+Dyl+Dyu+1,1)*sigma_int_z;
    
    epsr_x = ones(I+Dxl+Dxu+1,1)*epsr_z;
    epsr_y = ones(J+Dyl+Dyu+1,1)*epsr_z;
    
%        fprintf(1,'\n');
 end

%epsr_z

eps = epsr*epso;
mu  = mur*muo;

%auxilliary variables for conducting materials
conductive_aux.rho_x = zeros(1,I+Dxl+Dxu+1);
conductive_aux.rho_y = zeros(1,J+Dyl+Dyu+1);
conductive_aux.rho_z = zeros(1,K+Dzl+Dzu+1);

%Initialise the kappa matrices
kappa_x.E = ones(I+Dxl+Dxu+1,K+Dzl+Dzu+1);
kappa_x.H = ones(I+Dxl+Dxu+1,K+Dzl+Dzu+1);
kappa_y.E = ones(J+Dyl+Dyu+1,K+Dzl+Dzu+1);
kappa_y.H = ones(J+Dyl+Dyu+1,K+Dzl+Dzu+1);
kappa_z.E = ones(1,K+Dzl+Dzu+1);
kappa_z.H = ones(1,K+Dzl+Dzu+1);


%set up some values for determining the PML properties
pml_width_xl = (Dxl + 0.5)*delta.x;
pml_width_xu = (Dxu + 0.5)*delta.x;
pml_width_yl = (Dyl + 0.5)*delta.y;
pml_width_yu = (Dyu + 0.5)*delta.y;
pml_width_zl = (Dzl + 0.5)*delta.z;
pml_width_zu = (Dzu + 0.5)*delta.z;

sigma_max_xl = epso*c./sqrt(epsr_z*mur)*(n+1)*log(1/R0)/(pml_width_xl*2);
sigma_max_xu = epso*c./sqrt(epsr_z*mur)*(n+1)*log(1/R0)/(pml_width_xu*2);
sigma_max_yl = epso*c./sqrt(epsr_z*mur)*(n+1)*log(1/R0)/(pml_width_yl*2);
sigma_max_yu = epso*c./sqrt(epsr_z*mur)*(n+1)*log(1/R0)/(pml_width_yu*2);
%why was this Dxl+1 and Dxl+K?
sigma_max_zl = epso*c/sqrt(epsr_z(Dzl+1)*mur)*(n+1)*log(1/R0)/(pml_width_zl*2);
sigma_max_zu = epso*c/sqrt(epsr_z(max([Dzl+K,1]))*mur)*(n+1)*log(1/R0)/(pml_width_zu*2);
%sigma_max_zu = (c/sqrt(epsr_upperpml*mur_upperpml))*(n+1)*log(1/R0)/(pml_width_zu*2);

Nx = I+Dxl+Dxu;
Ny = J+Dyl+Dyu;
Nz = K+Dzl+Dzu;

%now setup the sigma_x
pmlvec_xl = (0:Dxl);
pmlvec_xu = (0:Dxu);
pmlvec_yl = (0:Dyl);
pmlvec_yu = (0:Dyu);
pmlvec_zl = (0:Dzl);
pmlvec_zu = (0:Dzu);

pmlvec_Se_xl = ((pmlvec_xl.' + 0.5)*delta.x/pml_width_xl).^n*sigma_max_xl;
pmlvec_Se_xu = ((pmlvec_xu.' + 0.5)*delta.x/pml_width_xu).^n*sigma_max_xu;
pmlvec_Se_yl = ((pmlvec_yl.' + 0.5)*delta.y/pml_width_yl).^n*sigma_max_yl;
pmlvec_Se_yu = ((pmlvec_yu.' + 0.5)*delta.y/pml_width_yu).^n*sigma_max_yu;
pmlvec_Se_zl = ((pmlvec_zl + 0.5)*delta.z/pml_width_zl).^n*sigma_max_zl;
pmlvec_Se_zu = ((pmlvec_zu + 0.5)*delta.z/pml_width_zu).^n*sigma_max_zu;
%save pmlvec_Se_xl pmlvec_Se_xl pmlvec_xl sigma_max_xl
%pmlvec_Se_zu = ((pmlvec_zu + 0.5)*delta.z/pml_width_zu).^n*eps*epsr_upperpml*sigma_max_zu;

pmlvec_Sh_xl = ((pmlvec_xl.' + 1.0)*delta.x/pml_width_xl).^n*sigma_max_xl;
pmlvec_Sh_xu = ((pmlvec_xu.' + 1.0)*delta.x/pml_width_xu).^n*sigma_max_xu;
pmlvec_Sh_yl = ((pmlvec_yl.' + 1.0)*delta.y/pml_width_yl).^n*sigma_max_yl;
pmlvec_Sh_yu = ((pmlvec_yu.' + 1.0)*delta.y/pml_width_yu).^n*sigma_max_yu;
pmlvec_Sh_zl = ((pmlvec_zl + 1.0)*delta.z/pml_width_zl).^n*sigma_max_zl;
pmlvec_Sh_zu = ((pmlvec_zu + 1.0)*delta.z/pml_width_zu).^n*sigma_max_zu;

%save pmlsig pmlvec_Se_xl pmlvec_Se_xu pmlvec_xu pml_width_xu
%sigma_max_xu pml_width_xl sigma_max_xl;

if Dxl>0
    kappa_x.E(1:(Dxl+1),:) = flipud(1 + ((pmlvec_xl.' + 0.5)*delta.x/pml_width_xl).^n*(kappa_max_z - 1));
end

if Dxu>0
    kappa_x.E((Nx - Dxu + 1):(Nx + 1),:) = 1 + ((pmlvec_xu.' + 0.5)*delta.x/pml_width_xu).^n*(kappa_max_z - 1);
end

if Dyl>0
    kappa_y.E(1:(Dyl+1),:) = flipud(1 + ((pmlvec_yl.' + 0.5)*delta.y/pml_width_yl).^n*(kappa_max_z - 1));
end

if Dyu>0
    kappa_y.E((Ny - Dyu + 1):(Ny + 1),:) = 1 + ((pmlvec_yu.' + 0.5)*delta.y/pml_width_yu).^n*(kappa_max_z - 1);
end

if Dzl>0
    kappa_z.E(1:(Dzl+1)) = fliplr(1 + (kappa_max(1) - 1)*((pmlvec_zl + 0.5)*delta.z/pml_width_zl).^n);
end

if Dzu>0
    kappa_z.E((Nz - Dzu + 1):(Nz + 1)) = 1 + (kappa_max(end) - 1)*((pmlvec_zu + 0.5)*delta.z/pml_width_zu).^n;
end

if Dxl>0
    kappa_x.H(1:Dxl,:) = flipud(1 + ((pmlvec_xl(1:(end-1)).'+ 1.0)*delta.x/pml_width_xl).^n*(kappa_max_z - 1));
end

if Dxu>0
    kappa_x.H((Nx - Dxu + 1):(Nx + 1),:) = 1 + ((pmlvec_xu.'+ 1.0)*delta.x/pml_width_xu).^n*(kappa_max_z - 1);
end

if Dyl>0
    kappa_y.H(1:Dyl,:) = flipud(1 + ((pmlvec_yl(1:(end-1)).'+ 1.0)*delta.y/pml_width_yl).^n*(kappa_max_z - 1));
end

if Dyu>0
    kappa_y.H((Ny - Dyu + 1):(Ny + 1),:) = 1 + ((pmlvec_yu.'+ 1.0)*delta.y/pml_width_yu).^n*(kappa_max_z - 1);
end

if Dzl>0
    kappa_z.H(1:Dzl) = fliplr(1 + (kappa_max(1) - 1)*((pmlvec_zl(1:(end-1))+ 1.0)*delta.z/pml_width_zl).^n);
end

if Dzu>0
    kappa_z.H((Nz - Dzu + 1):(Nz + 1)) = 1 + (kappa_max(end) - 1)*((pmlvec_zu+ 1.0)*delta.z/pml_width_zu).^n;
end

save kappa kappa_x kappa_y kappa_z kappa_max_z;

%now populate the sigma_x matrix
Dx = Dxl;
pmlvec_Se = pmlvec_Se_xl;
pmlvec_Sh = pmlvec_Sh_xl;
if Dx ~= 0
    for i=1:(Dx+1)
	sigma_x.E(i,:) = pmlvec_Se(Dx+1-i+1,:);    
	if i < Dx + 1
	    sigma_x.H(i,:) = pmlvec_Sh(Dx+1-i,:);    
	end
    end
end
    
Dx = Dxu;
pmlvec_Se = pmlvec_Se_xu;
pmlvec_Sh = pmlvec_Sh_xu;
if Dx ~= 0
    for i=(Nx - Dx + 1):(Nx + 1)
        sigma_x.E(i,:) = pmlvec_Se(i-Nx+Dx,:);
        sigma_x.H(i,:) = pmlvec_Sh(i-Nx+Dx,:);
    end
end

%now setup the sigma_y
%pmlvec = delta*(0:Dy);
%pmlvec_Se = ((pmlvec + 0.5)*delta/pml_width_y).^n*eps*sigma_max_y;
%pmlvec_Sh = ((pmlvec + 0.5)*delta/pml_width_y).^n*mu*sigma_max_y;
Dy = Dyl;
pmlvec_Se = pmlvec_Se_yl;
pmlvec_Sh = pmlvec_Sh_yl;

%now populate the sigma_y matrix
if Dy ~= 0
    for j=1:(Dy+1)
	sigma_y.E(j,:) = pmlvec_Se(Dy+1-j+1,:);    
	if j < Dy + 1
	    sigma_y.H(j,:) = pmlvec_Sh(Dy+1-j,:);    
	end
    end
end

Dy = Dyu;
pmlvec_Se = pmlvec_Se_yu;
pmlvec_Sh = pmlvec_Sh_yu;

if Dy ~= 0
    for j=(Ny - Dy + 1):(Ny + 1)
        sigma_y.E(j,:) = pmlvec_Se(j-Ny+Dy,:);
        sigma_y.H(j,:) = pmlvec_Sh(j-Ny+Dy,:);
    end
end

%now setup the sigma_z
%pmlvec = delta*(0:Dz);
%pmlvec_Se = ((pmlvec + 0.5)*delta/pml_width_y).^n*eps*sigma_max_y;
%pmlvec_Sh = ((pmlvec + 0.5)*delta/pml_width_y).^n*mu*sigma_max_y;

Dz = Dzl;
pmlvec_Se = pmlvec_Se_zl;
pmlvec_Sh = pmlvec_Sh_zl;

%now populate the sigma_z matrix
if Dz ~= 0
    for k=1:(Dz+1)
	sigma_z.E(k) = pmlvec_Se(Dz+1-k+1);    
	if k < Dz + 1
	    sigma_z.H(k) = pmlvec_Sh(Dz+1-k);    
	end
    end
end


Dz = Dzu;
pmlvec_Se = pmlvec_Se_zu;
pmlvec_Sh = pmlvec_Sh_zu;

if Dz ~= 0
    for k=(Nz - Dz + 1):(Nz + 1)
        sigma_z.E(k) = pmlvec_Se(k-Nz+Dz);
        sigma_z.H(k) = pmlvec_Sh(k-Nz+Dz);
    end
end


%dispersive parameters
dispersive_aux.alpha = 4./(vc_vec*dt + 2);
dispersive_aux.beta = (vc_vec*dt - 2)./(vc_vec*dt + 2);
dispersive_aux.gamma = 2*epso*wp_vec.^2*dt^2./(vc_vec*dt + 2);
dispersive_aux.kappa_x = kappa_x.E;
dispersive_aux.kappa_y = kappa_y.E;
dispersive_aux.kappa_z = kappa_z.E;
dispersive_aux.sigma_x = sigma_x.E;%(:,1);
dispersive_aux.sigma_y = sigma_y.E;%(:,1);
dispersive_aux.sigma_z = sigma_z.E;


%now setup the actual C and D parameters
%save sigma sigma_z sigma_y sigma_x;
%initialise the param matrices
Cax =  zeros(I+Dxl+Dxu+1,K+Dzl+Dzu+1);
Cay =  zeros(J+Dyl+Dyu+1,K+Dzl+Dzu+1);
Caz =  zeros(1,K+Dzl+Dzu+1);

Dax =  zeros(I+Dxl+Dxu+1,K+Dzl+Dzu+1);
Day =  zeros(J+Dyl+Dyu+1,K+Dzl+Dzu+1);
Daz =  zeros(1,K+Dzl+Dzu+1);

Cbx =  zeros(I+Dxl+Dxu+1,K+Dzl+Dzu+1);
Cby =  zeros(J+Dyl+Dyu+1,K+Dzl+Dzu+1);
Cbz =  zeros(1,K+Dzl+Dzu+1);

Dbx =  zeros(I+Dxl+Dxu+1,K+Dzl+Dzu+1);
Dby =  zeros(J+Dyl+Dyu+1,K+Dzl+Dzu+1);
Dbz =  zeros(1,K+Dzl+Dzu+1);


%inside the free space region
%sigma_int = 1e4;

ax = dt*(kappa_x.E.*sigma_int_x + epsr_x.*sigma_x.E)./(2*epsr_x*epso.*kappa_x.E) + dt*dt*sigma_int_x.*sigma_x.E./(4*epso*eps*kappa_x.E);
bx = (ones(I+Dxl+Dxu+1,1)*dispersive_aux.gamma)./(4*epsr_x*epso);
cx = dt*sigma_x.E.*(ones(Dxl+Dxu+I+1,1)*dispersive_aux.gamma)./(2*kappa_x.E*epso^2.*epsr_x);
Cax = (1 - ax - cx)./(1 + ax + bx);
Cbx = dt./(epsr_x*epso.*kappa_x.E*delta.x)./(1 + ax + bx);
Ccx = (ones(Dxl+Dxu+I+1,1)*dispersive_aux.gamma)./(4*epso*epsr_x)./(1 + ax + bx);



ax = sigma_x.H.*dt./(2*epso*kappa_x.H);
Dax = (1 - ax)./(1 + ax);
Dbx = dt./(muo*kappa_x.H*delta.x)./(1 + ax);

ay = dt*(kappa_y.E.*sigma_int_y + epsr_y.*sigma_y.E)./(2*epsr_y*epso.*kappa_y.E) + dt*dt*sigma_int_y.*sigma_y.E./(4*epso*epso*epsr_y.*kappa_y.E);
by = (ones(J+Dyl+Dyu+1,1)*dispersive_aux.gamma)./(4*eps);
cy = dt*sigma_y.E.*(ones(Dyl+Dyu+J+1,1)*dispersive_aux.gamma)./(2*kappa_y.E*epso^2.*epsr_y);
Cay = (1 - ay - cy)./(1 + ay + by);
Cby = dt./(epsr_y*epso.*kappa_y.E*delta.y)./(1 + ay + by);
Ccy = (ones(Dyl+Dyu+J+1,1)*dispersive_aux.gamma)./(4*epso*epsr_y)./(1 + ay + by);

ay = sigma_y.H*dt./(2*epso*kappa_y.H);
Day = (1 - ay)./(1 + ay);
Dby = dt./(muo*kappa_y.H*delta.y)./(1 + ay);

az = dt*(kappa_z.E.*sigma_int_z + epsr_z.*sigma_z.E)./(2*epsr_z*epso.*kappa_z.E) + dt*dt*sigma_int_z.*sigma_z.E./(4*epso*epsr_z*epso.*kappa_z.E);
bz = dispersive_aux.gamma./(4*epsr_z*epso);
cz = dt*sigma_z.E.*dispersive_aux.gamma./(2*kappa_z.E*epso^2.*epsr_z);
Caz = (1 - az - cz)./(1 + az + bz);
Cbz = dt./(epsr_z*epso.*kappa_z.E*delta.z)./(1 + az + bz);
Ccz = dispersive_aux.gamma./(4*epso*epsr_z)./(1 + az + bz);

save Cvars Caz Cbz Ccz az bz cz;

az = sigma_z.H*dt./(2*epso*kappa_z.H);
Daz = (1 - az)./(1 + az);
Dbz = dt./(muo*kappa_z.H*delta.z)./(1 + az);

conductive_aux.rho_x = dt*sigma_int_x.*sigma_x.E/(2*epso);
conductive_aux.rho_y = dt*sigma_int_y.*sigma_y.E/(2*epso);
conductive_aux.rho_z = dt*sigma_int_z.*sigma_z.E/(2*epso);

sigma.sigma_x = sigma_x;
sigma.sigma_y = sigma_y;
sigma.sigma_z = sigma_z;

if isempty(multilayer)
    Cax = Cax(:,1).';
    Cay = Cay(:,1).';
    Cbx = Cbx(:,1).';
    Cby = Cby(:,1).';
    Ccx = Ccx(:,1).';
    Ccy = Ccy(:,1).';
    Dax = Dax(:,1).';
    Day = Day(:,1).';
    Dbx = Dbx(:,1).';
    Dby = Dby(:,1).';
    
    conductive_aux.rho_x = conductive_aux.rho_x(:,1).';
    conductive_aux.rho_y = conductive_aux.rho_y(:,1).';
end

C.Cax = Cax;
C.Cay = Cay;
C.Caz = Caz;
C.Cbx = Cbx;
C.Cby = Cby;
C.Cbz = Cbz;
C.Ccx = Ccx;
C.Ccy = Ccy;
C.Ccz = Ccz;

D.Dax = Dax;
D.Day = Day;
D.Daz = Daz;
D.Dbx = Dbx;
D.Dby = Dby;
D.Dbz = Dbz;

save dispersive_aux dispersive_aux C;

%now calculate the b variables for free space
freespace.Cby = dt/(eps*delta.y);
freespace.Cbz = dt/(eps*delta.z);
freespace.Cbx = dt/(eps*delta.x);

freespace.Dby = dt/(mu*delta.y);
freespace.Dbz = dt/(mu*delta.z);
freespace.Dbx = dt/(mu*delta.x);



