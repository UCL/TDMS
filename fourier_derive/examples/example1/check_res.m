dat=load('out.txt');

%N = size(dat,1)/2;
N = 32;

xc = dat(1:N,1);
xs = dat((N+1):(N+1+N-1),1);

g = zeros(2,2,N);
for i=1:2
    for j=1:2
	g(i,j,:) = dat(1:N,2);
    end
end
%g_d = FFT_derivative_3d(g,1/N,0.5/N,3);
g_d = numerical_derivative(dat(1:N,2),0.5);
figure(1);clf;
%plot(xc,cos(2*pi*xc)+sin(2*pi*2*xc));
%hold on;
%plot(xc,dat(1:N,2),'r--');
plot(xs,-2*pi*sin(2*pi*xs)+4*pi*cos(2*pi*2*xs));
hold on;
plot(xs,dat((N+1):(N+1+N-1),2),'r--');
plot(xs,squeeze(g_d),'k.');

figure(2);clf;
plot(xc,squeeze(g(1,1,:)));

dk1 =  dat((2*N+1):(2*N+1+N-1),2) + sqrt(-1)*dat((2*N+1):(2*N+1+N-1),3);
dat=load('G_array');
dk2 = dat.G_array.';


figure(3);clf;
plot(real(dk1));
hold on;
plot(real(dk2),'r');
plot(imag(dk1),'--');
hold on;
plot(imag(dk2),'r--');
