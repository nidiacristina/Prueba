clear all 
clc

N=10;
M=10;

xa=0;
ya=0;
xb=4;
yb=-1;
xc=7;
yc=3;
xd=1;
yd=4;


f1= @(x) x*(yb-ya)/(xb-xa);
f2= @(x) x*(xc-xb)/(yc-yb);
f3= @(x) x*(yd-yc)/(xd-xc);
f4= @(x) x*(xa-xd)/(ya-yd);

f1_prime= @(x) (yb-ya)/(xb-xa);
f2_prime= @(x) (xc-xb)/(yc-yb);
f3_prime= @(x) (yd-yc)/(xd-xc);
f4_prime= @(x) (xa-xd)/(ya-yd);


[X_el,Y_el]=elliptic_mapping(N,M,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd,0);


figure
mesh(X_el,Y_el,zeros(M,N))
view(2);
axis('equal');
