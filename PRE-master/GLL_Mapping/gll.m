function [x,w] = gll(N)
% Calculate the N GLL points and interpolation weights
maxIter=20;

x=-1.0*cos(pi*(0:N-1)/(N-1))';

for ii=1:maxIter
   xold=x; 
   
   Lk=ones(N,1);
   Lkp1=xold;
   
   for kk=2:N-1
      Lkm1=Lk;
      Lk=Lkp1;
      Lkp1=((2*kk-1)*xold.*Lk - (kk-1)*Lkm1)/kk;
   end
   
   x = xold - (xold.*Lkp1-Lk)./(N*Lkp1);
   
   if(max(abs(x-xold)) < 2*eps) 
       break;
   end
end

%display(num2str(ii));

if(ii == maxIter)
   display('MAXIMUM ITERATION REACHED IN GLL!!'); 
end

w= 2.0./((N-1)*N*Lkp1.^2);




