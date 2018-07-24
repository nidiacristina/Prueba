clear all 
clc

N=25;

f1= @(x) sin(2*x)/4;
f2= @(x) sin(2*x)/4;
f3= @(x) sin(2*x)/4;
f4= @(x) sin(2*x)/4;

f1_prime= @(x) cos(2*x)/2;
f2_prime= @(x) cos(2*x)/2;
f3_prime= @(x) cos(2*x)/2;
f4_prime= @(x) cos(2*x)/2;

xa=0;
ya=0;
xb=2*pi;
yb=0;
xc=2*pi;
yc=2*pi;
xd=0;
yd=2*pi;

% f1= @(x) 0;
% f2= @(x) 0;
% f3= @(x) 0;
% f4= @(x) 0;
% 
% f1_prime= @(x) 0;
% f2_prime= @(x) 0;
% f3_prime= @(x) 0;
% f4_prime= @(x) 0;
% 
% xa=0;
% ya=0;
% xb=1;
% yb=0;
% xc=1;
% yc=1;
% xd=0;
% yd=1;

% [X,w]=gll_arclengthmethod(N,f1,f1_prime,yc,xc,yb,xb);
[X_GH,Y_GH]=gordon_hall_mapping(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd);
%[X,w] = gll(N);
figure
mesh(X_GH,Y_GH,zeros(N,N))
view(2);
axis('equal');