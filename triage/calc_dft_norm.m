function [E_norm] = calc_dft_norm(inputfile)

%addpath('/home/pmunro/code/FDTD/dispersive1.2_dg2d_devel/');
[epso , muo , c] = import_constants;
%rmpath('/home/pmunro/code/FDTD/dispersive1.2_dg2d_devel/');


[fid_input,message] = fopen(inputfile,'r');

%check if file was opened successfully
if fid_input== -1
    error(sprintf('File %s could not be opened for reading',input_file));
end

%proceed to_l read in config information
current_line = fgets(fid_input);
while current_line ~= -1
    if ~isempty(findstr('material_file',current_line)) 
	if isempty(material_file) | ((~isempty(material_file)) & (findstr('material_file',current_line)>findstr('%',current_line)))
	    eval(current_line);
	else 
	    fprintf(1,'material_file specified as argument and in %s, using %s\n',input_file,material_file);
	end
    else
	eval(current_line);
    end
    
    current_line = fgets(fid_input);
end
fclose(fid_input);
%input_file_tfoc_master_compare_50;


if 0
    wavelengthwidth = 100e-9;
    dat_lambda = load('/home/pmunro/research/colleagues/philip_wijesinghe/temporal_microsc/code/PSTD/11/compare_tf_cw/lambda_vec_sample_50');
    lambda = 800e-9;
    k_vec = 2*pi./dat_lambda.lambda_vec_sample;
    f_ex_vec=2.997924580105029e+08*k_vec/2/pi;
    dt = 2/sqrt(2)/pi*(lambda/6)/(3e8/1.0)*.95;
    start_tind = 0;
    Nt = 211500;
    epsr = [1^2];
end

start_tind = 0;
refractive_index = sqrt(real(epsr(1)));
f_an =  asin( 2*pi/lambda*2.997924580105029e+08*dt/2)/(pi*dt);
omega_an = 2*pi*f_an;
lambda_an = c/(f_an*refractive_index);
hwhm = lambda_an^2/((c/refractive_index)*wavelengthwidth)*2*sqrt(log(2)/pi);
to_l = hwhm*sqrt(log(1e8)/pi);
E_norm = zeros(size(f_ex_vec));
f_max=max(f_ex_vec);
Np=floor(1./(2.5*dt*f_max));

dtp = Np*dt;

Npe = 0;
for tind = start_tind:(Nt-1)
    if( mod(tind-start_tind,Np) == 0)
      Npe=Npe+1;
    end
end
Nsteps = Npe;
fte_acc = zeros(1,numel(0:(Nt-1)));

for tind=0:(Nt-1)
    time_E = (tind + 1)*dt;
    fte = real((-1.*sqrt(-1))*exp(-sqrt(-1)*(omega_an*(time_E - to_l))))*exp(-1.*pi*((time_E - to_l)/(hwhm)).^2);
    fte_acc(tind+1) = fte;
    
    if( mod(tind-start_tind,Np) == 0)
	for ifx=1:numel(f_ex_vec)
	    omega = f_ex_vec(ifx)*2*pi;
	    E_norm(ifx) = E_norm(ifx) + fte*exp( (omega*( (tind+1))*dt) * sqrt(-1)) * 1./( Npe);
	end
    end
end


