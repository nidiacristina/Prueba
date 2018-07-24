clear all 
clc

N=40;
M=10;

f1= @(x) sin(2*(x-3))/4;
f2= @(x) sin(2*x)/4;
f3= @(x) 0;
f4= @(x) 0;

f1_prime= @(x) cos(2*(x-3))/2;
f2_prime= @(x) cos(2*x)/2;
f3_prime= @(x) 0;
f4_prime= @(x) 0;

xa=3;
ya=0;
xb=2*pi;
yb=0;
xc=2*pi;
yc=2*pi;
xd=3;
yd=2*pi;
% [X,w]=gll_arclengthmethod(N,f1,f1_prime,yc,xc,yb,xb);
%[X_GH,Y_GH]=gordon_hall_mapping(N,M,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
%f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd);
[X_el,Y_el]=elliptic_mapping(N,M,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd,0);
% %[X,w] = gll(N);
% [s,t]=control_mapping(N, 'BOUNDARY_ORTHOGONALITY');
figure
%plot(X_el,Y_el,'*')
% mesh(X_el,Y_el,zeros(N,N))
mesh(X_el,Y_el,zeros(M,N))
view(2);
axis('equal');

% plot(s,t,'*')