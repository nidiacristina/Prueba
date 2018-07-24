clear all 
clc

N=30;
M=15;
n_iter=1;

filename='Laplace_sin_25';

f1= @(x) sin(x)/4;
% f1= @(x) 0;
f2= @(x) sin(x)/2;
f3= @(x) 0;
f4= @(x) 0;

f1_prime= @(x) cos(x)/4;
% f1_prime= @(x) 0;
f2_prime= @(x) cos(x)/2;
f3_prime= @(x) 0;
f4_prime= @(x) 0;

xa=0;
ya=0;
xb=2*pi;
yb=0;
xc=4*pi;
yc=0;
xd=4*pi;
yd=pi;
xe=2*pi;
ye=pi;
xf=0;
yf=pi;
xg=0;
yg=2*pi;
xh=2*pi;
yh=2*pi;
xi=4*pi;
yi=2*pi;

% [X_GH,Y_GH]=gordon_hall_mapping_2(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
% f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd);

[X_el_1,Y_el_1]=elliptic_mapping(N,M,f3,f3_prime,f1,f1_prime,f1,f1_prime,...
f3,f3_prime,xa,ya,xb,yb,xe,ye,xf,yf,n_iter);
[X_el_2,Y_el_2]=elliptic_mapping(N,M,f3,f3_prime,f3,f3_prime,f1,f1_prime,...
f1,f1_prime,xb,yb,xc,yc,xd,yd,xe,ye,n_iter);
[X_el_3,Y_el_3]=elliptic_mapping(N,M,f1,f1_prime,f1,f1_prime,f3,f3_prime,...
f3,f3_prime,xf,yf,xe,ye,xh,yh,xg,yg,n_iter);
[X_el_4,Y_el_4]=elliptic_mapping(N,M,f1,f1_prime,f3,f3_prime,f3,f3_prime,...
f1,f1_prime,xe,ye,xd,yd,xi,yi,xh,yh,n_iter);


X_el=[X_el_3(1:end,1:end-1) X_el_4; X_el_1(2:end,1:end-1) X_el_2(2:end,1:end)];
Y_el=[Y_el_3(1:end,1:end-1) Y_el_4; Y_el_1(2:end,1:end-1) Y_el_2(2:end,1:end)];

figure
%mesh(X_el,Y_el,zeros(2*M-1,2*N-1))
mesh(X_el_1,Y_el_1,zeros(M,N))
view(2);
axis('equal');

%save(filename,'X_el','X_el_1','X_el_2','X_el_3','X_el_4','Y_el','Y_el_1','Y_el_2','Y_el_3','Y_el_4')