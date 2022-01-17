lambda = 405e-9;

z = linspace(0,8*lambda,100);
z2 = linspace(0,8*lambda,8*2.5);
dz2 = z2(2)-z2(1);

dz=z(2)-z(1);
z = [-dz z];


k=2*pi/lambda;
W=20*lambda;
Hy = sin(k*(z-2*dz)).*exp( -((z-2*lambda)/W).^2 );

dz_s = ((z(1:(end-1))))+dz/2;
dz_c = ((z(1:(end-1))));

dhy_s = k*cos(k*(dz_s-2*dz)).*exp( -((dz_s-2*lambda)/W).^2 )-2/W*((dz_s-2*lambda)/W).*sin(k*(dz_s-2*dz)).*exp( -((dz_s-2*lambda)/W).^2 );
dhy_s(find((dz_s/1e-6-0.5237)<=0))=0;

dhy_c = k*cos(k*(dz_c-2*dz)).*exp( -((dz_c-2*lambda)/W).^2 )-2/W*((dz_c-2*lambda)/W).*sin(k*(dz_c-2*dz)).*exp( -((dz_c-2*lambda)/W).^2 );


in = 18;
dz_c = [dz_c(1:in) dz_c(in) dz_c( (in+1):end)];
dhy_c =[dhy_c(1:in) dhy_c(in) dhy_c( (in+1):end)];
dhy_c(1:in)=0;
%dhy_c(find((dz_c/1e-6-0.5237)<=0))=0;

Hy(1:18)=0;

Hy2 = sin(k*(z2-2*dz)).*exp( -((z2-2*lambda)/W).^2 );
Hy2(1:4) = 0;

dHy_1 = diff(Hy)/dz;

Hyt = zeros(2,2,numel(Hy));
for i=1:2
    for j=1:2
	Hyt(i,j,:)=Hy;
    end
end
dHy_2 = FFT_derivative_3d(Hyt,dz,1/2*dz,3);
dHy_2 = squeeze(dHy_2(1,1,1:(end-1)));

Hy2t = zeros(2,2,numel(Hy2));
for i=1:2
    for j=1:2
	Hy2t(i,j,:)=Hy2;
    end
end
dHy2_2 = FFT_derivative_3d(Hy2t,dz2,1/2*dz2,3);
dHy2_2 = squeeze(dHy2_2(1,1,1:(end-1)));


dHy_3 = FFT_derivative_3d(Hyt,dz,0,3);
dHy_3 = squeeze(dHy_3(1,1,1:(end-1)));

fs = 16;lw = 2;
figure(1);clf;

subplot(2,2,1);
h1 = plot(z/1e-6-0.5237,Hy,'.-');axis tight;grid on;
ht = text(-.4,.8,'a)');

set(gca,'FontSize',fs);
set(ht,'Units','Normalized','FontSize',fs,'Interpreter','latex');
apos = get(ht,'Position');

dleg = 0.;
ax=axis;ax(2)=1.5-0.5237+dleg;axis([ax(1:2) ax(3:4)*1.05]);
hy = ylabel('$H_y (a.u.)$');
hx = xlabel('$z (\mu m)$');
set([hy hx ],'Interpreter','latex');

subplot(2,2,4);
h3 = plot((z(1:(end-1))+dz/2)/1e-6-0.5237,dHy_2/1e6,'.-');axis tight;grid on;
set(gca,'FontSize',fs);
ax=axis;ax(2)=1.5-0.5237+dleg;axis([ax(1:2) ax(3:4)*1.05]);
ax2=axis;hold on;
ha4 = plot(dz_s/1e-6-0.5237,dhy_s/1e6,'rx');

hy = ylabel('$\partial H_y/\partial z (a.u./\mu m)$');
hx = xlabel('$z (\mu m)$');
ht = title('    Fourier, staggered');
set([hy hx ht],'Interpreter','latex');

ht = text(-.4,.8,'d)');
set(ht,'Units','Normalized','FontSize',fs,'Interpreter','latex');
set(ht,'Position',apos);

set(gca,'XTick',[-0.5 dz/2*1e6 0.5]);
set(gca,'XTickLabel',{'-0.5','','0.5'});

htdelta = text(dz/2*1e6-.2,-19.5,'$\Delta_{xz}/2$','Interpreter','latex','FontSize',fs);


subplot(2,2,3);
h4 = plot((z(1:(end-1)))/1e-6-0.5237,dHy_3/1e6,'.-');axis tight;grid on;
ax=axis;ax(2)=1.5-0.5237+dleg;axis([ax(1:2) ax(3:4)*1.05]);set(gca,'FontSize',fs);hold on;
ha3 = plot(dz_c/1e-6-0.5237,dhy_c/1e6,'rx');
hy = ylabel('$\partial H_y/\partial z (a.u./\mu m)$');
hx = xlabel('$z (\mu m)$');

ht = title('    Fourier, collocated');
set([hy hx ht],'Interpreter','latex');
ht = text(-.4,.8,'c)');
set(ht,'Units','Normalized','FontSize',fs,'Interpreter','latex');
set(ht,'Position',apos);



subplot(2,2,2);
h2 = plot((z(1:(end-1))+dz/2)/1e-6-0.5237,dHy_1/1e6,'-');axis tight;grid on;
hold on;
h2m = plot((z(1:(end-1))+dz/2)/1e-6-0.5237,dHy_1/1e6,'.');axis tight;grid on;
ha1 = plot(dz_s/1e-6-0.5237,dhy_s/1e6,'rx');

set(gca,'FontSize',fs);
ax=axis;ax(2)=1.5-0.5237+dleg;axis([ax(1:2) ax(3:4)*1.05]);
axis(ax2);

set(gca,'XTick',[-0.5 dz/2*1e6 0.5]);
set(gca,'XTickLabel',{'-0.5','','0.5'});

htdelta = text(dz/2*1e6-.2,-19.5,'$\Delta_{xz}/2$','Interpreter','latex','FontSize',fs);

hy = ylabel('$\partial H_y/\partial z (a.u./\mu m)$');
hx = xlabel('$z (\mu m)$');
ht = title('    Central differences');
set([hy hx ht],'Interpreter','latex');


set([h1 h2 h3 h4],'LineWidth',lw);
set([h1 h2m h3 h4],'MarkerSize',15);
set([ ha1 ha3 ha4],'MarkerSize',12);
set(gca,'FontSize',fs);

hl = legend([h2m ha1],{'Numeric','Analytic'},'Location','NorthEast');
pos = get(hl,'Position');
set(hl,'Position',pos + [0.06 .03 0 0]);

set(hl,'Interpreter','latex');set(hl,'FontSize',12)
ht = text(-.4,.8,'b)');
set(ht,'Units','Normalized','FontSize',fs,'Interpreter','latex');
set(ht,'Position',apos);

%print -depsc2 example_derivatives.eps;


figure(2);clf;
subplot(2,1,1);
plot(z,Hy);hold on;
plot(z2,Hy2,'r.');
