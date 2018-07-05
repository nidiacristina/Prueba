function [y_ap] = Runge_Kutta_4(f,x0,xf,y0,h)

%f = @(x,y) 2*y*(1-(y/100)) %f(x,y) - EDO
%ana=@(x,y) -2/4*x.^4+4*x.^3-10*x.^2+8.5*x+1; %Analítica
%x0=0; %Condiiones iniciales
%y0=1; %Condiiones iniciales
%xo=0;
%xf=10;

%xo=x0;
%xf=0.5;
%h=0.25;

lx=[x0:h:xf];
y_ap=zeros(size(lx));

%Calcular aproximación numérica
x=x0;
y=y0;
for i=1:length(lx)
    k1=f(x,y);
    k2=f(x+0.5*h,y+0.5*h*k1);
    k3=f(x+0.5*h,y+0.5*k2*h);
    k4=f(x+h,y+k3*h);
    y1=y+1/6*(k1+2*k2+2*k3+k4)*h;
    y_ap(i)=y1;
    
    x=x+h;
    y=y1;
end

y_ap=[y0 y_ap(1:(end-1))];
%Clacular analítica
%lx2=[xo:h/10:xf];
%y_ana=ana(lx2,0);

%{
%Graficar
y_ana2=ana(lx,0);
error=abs(y_ana2-y_ap)./y_ana2;

subplot(3,1,[1:2])
plot(lx2,y_ana,'b')
hold on
plot(lx,y_ap,'o-r')
%legend('Analitica','Aproximacion')
ylabel('y')
title('RK4')
subplot(3,1,3)
%plot(lx,error,'k')
semilogy(lx,error,'k')
ylabel('Error')
xlabel('x')
max_er=max(error);
title(['Max error-->',num2str(max_er)])
%}
end
