%Prueba matriz de diferenciaci[on
clc
clear 
close all

zlims=[0 2*pi];
n=3;
mz=500;
%Matriz dif finitas


%AD=global_spectral_differentiation_matrix( zlims, n, mz);

x=linspace(zlims(1),zlims(2),n*mz)';
dx=x(2)-x(1);
AD=full(gallery('tridiag',n*mz,-1,1,0))./(dx);
AD(1,2)=0;
AD(1,1)=1;
ana=sin(x);

%Vector mano derecha
b=cos(x);

num=AD*ana;

%PLot
plot(x,b)
hold on
plot(x,num)



