datt=load('out/fdtd_cyl_tdms');
datp=load('out/fdtd_cyl');
%datt=load('finalgrid_01');
%datp=load('finalgrid');
%datt=datt.fdtdgrid;
%datp=datp.fdtdgrid;

fields=fieldnames(datt);

err = zeros(1,numel(fields));
for i=1:numel(fields)
    if ~isempty(getfield(datt,fields{i})) & ~isempty(getfield(datp,fields{i}))
	A = getfield(datt,fields{i});
	B = getfield(datp,fields{i});
	if sqrt(sum(abs(A(:)).^2))>0
	    err(i) = sqrt(sum(abs(A(:)-B(:)).^2))/sqrt(sum(abs(A(:)).^2));
	else
	    err(i) = sqrt(sum(abs(A(:)-B(:)).^2));
	end
    end
end

%Ex = datp.Exy+datp.Exz;
%Ey = datp.Eyx+datp.Eyz;
%Ez = datp.Ezy+datp.Ezx;
%Hx = datp.Hxy+datp.Hxz;
%Hy = datp.Hyx+datp.Hyz;
%Hz = datp.Hzy+datp.Hzx;

%maxfield=max([max(abs(Ex(:))) max(abs(Ey(:))) max(abs(Ez(:))) max(abs(Hx(:))) max(abs(Hy(:))) max(abs(Hz(:)))])
