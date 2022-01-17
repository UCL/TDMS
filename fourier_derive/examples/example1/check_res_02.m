N = 33;
delta = 0.5;

x = linspace(0,1,N+1);xhr = linspace(0,1,500);xhr = xhr(1:(end-1));
x = x(1:(end-1));
dx = diff(x(1:2));



%fvec  = round(rand(1,Ns)*10);
%fvec = [1 2 4 8 7];% 2 4 8 16];
%fvec = rand(1,200)*12;
fvec = 14;
Ns = numel(fvec);
avec = rand(1,Ns);
pvec = rand(1,Ns)*2*pi;


g_a = zeros(size(x));
g_d_a = zeros(size(x));g_d_a_hr = zeros(size(xhr));g_a_hr = zeros(size(xhr));
W = .15;
xg = 0.5;
mask = double(x>.4);maskhr = double(xhr>.4);
for i=1:Ns
    g_a = g_a + exp(-(x-xg).^2/W^2).*avec(i).*sin(2*pi*fvec(i)*(x-pvec(i)));
    g_a_hr = g_a_hr + exp(-(xhr-xg).^2/W^2).*avec(i).*sin(2*pi*fvec(i)*(xhr-pvec(i)));
    g_d_a = g_d_a + exp(-(x-xg).^2/W^2).*avec(i).*(2*pi*fvec(i)).*cos(2*pi*fvec(i)*(x-pvec(i)));
    g_d_a_hr = g_d_a_hr + exp(-(xhr-xg).^2/W^2).*avec(i).*(2*pi*fvec(i)).*cos(2*pi*fvec(i)*(xhr-pvec(i)));
end
g_d_a = g_d_a -2*(x-xg)/W^2.*g_a;
g_d_a_hr = g_d_a_hr -2*(xhr-xg)/W^2.*g_a_hr;

g_a = g_a.*mask;
g_d_a = g_d_a.*mask;
g_a_hr = g_a_hr.*maskhr;
g_d_a_hr = g_d_a_hr.*maskhr;




G3_a = zeros(2,2,numel(g_a));
for i=1:2
    for j=1:2
	G3_a(i,j,:) = g_a;
    end
end
G3_d_a = FFT_derivative_3d(G3_a,1/N,delta/N,3);


g_d_n =  numerical_derivative(g_a,delta);

figure(1);clf;
plot(x,g_d_a,'*');hold on;
plot(xhr,g_d_a_hr,'-');

plot(x+delta*dx,g_d_n,'r.');
plot(x+delta*dx,squeeze(G3_d_a(1,1,:)),'kd');

