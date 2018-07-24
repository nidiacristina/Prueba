
%Return the barycentric interpolation weights and differentiation
%matrices for a given grid, x, using Berrut & Trefethen SIAM Rev. 2004
function [omega,D1,D2]=baryc(x)
N=length(x);

omega=zeros(N,1);
D1=zeros(N,N);
D2=zeros(N,N);

%Compute the barycentric interpolation weights
for j=1:N
    
    xj=x(j);
    eta=omega(j);
    for k=[1:j-1,j+1:N]
        eta=eta+log(abs(xj-x(k)));
    end
    omega(j)=((-1)^(j))*exp(-eta);

end

%Compute the 1st derivative matrix
for j=1:N

    xj=x(j);
    for k=[1:j-1,j+1:N]
       D1(j,k)=(omega(k)/omega(j))/(xj-x(k)); 
    end
    
    D1(j,j)=-1.0*sum(D1(j,:));
end

% for j=1:N
%     D1(j,j)=-1.0*sum(D1(j,:));
% endD-D1
    
%Compute the 2nd derivative matrix
for j=1:N
    
    xj=x(j);
    for k=[1:j-1,j+1:N]
        D2(j,k)=-2.0*D1(j,k)*(1.0/(xj-x(k))-D1(j,j));
    end
    
    D2(j,j)=-1.0*sum(D2(j,:));
end
