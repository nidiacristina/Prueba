function [X,w]=gll_arclengthmethod_horizontal(N,f,f_prime,xa,ya,xb,yb)
% Calculate the N GLL points and interpolation weights using arc length
% method

[x,w]=gll(N);
x=1.+x; %GLL points for segment [0,2]
L=arclength(f_prime,xa,xb); %total arc length

x_L=x.*L/2; %GLL points for segment [0,L] 

X=zeros(N,2); %Coordinates of GLL points 
X(1,1)=xa;
X(1,2)=ya;
X(N,1)=xb;
X(N,2)=yb;
arclen= @(z) arclength(f_prime,xa,z);

for k=2:N-1
    arcf= @(y) arclen(y)-x_L(k);
    X(k,1)=fzero(arcf, X(k-1,1));
    X(k,2)=f(X(k,1)-xa)+ya;
end
    
end

