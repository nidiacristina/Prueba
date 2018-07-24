clear all 
clc

N_1=24;
M_1=8;
N_2=15;
M_2=15;



f1= @(x) 8*sin(x/10);
% f1= @(x) 0;
f2= @(x) 0;
f3= @(x) 0;
f4= @(x) 0;

f1_prime= @(x) 8/10*cos(x/10);
% f1_prime= @(x) 0;
f2_prime= @(x) 0;
f3_prime= @(x) 0;
f4_prime= @(x) 0;

xa=0;
ya=0;
xb=100*pi;
yb=0;
xc=100*pi;
yc=100;
xd=0;
yd=100;

[X_grid,Y_grid]=general_grid_generation(N_1,M_1,N_2,M_2,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd,'equidistant','GH');

figure
mesh(X_grid,Y_grid,zeros((M_1-1)*M_2,(N_1-1)*M_2))
view(2);
axis('equal');