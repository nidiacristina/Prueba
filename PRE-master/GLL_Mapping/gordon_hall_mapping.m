function [X_GH,Y_GH]=gordon_hall_mapping(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd)
%Calculate the position of the points of the mapping with Gordon Hall
%algorithm
%                   f3
%       xd  x - - - - - - - x xc
%           |               |
%           |               |
%        f3 |               | f2
%           |               |
%           |               |
%           |               |
%       xa  x - - - - - - - x xb
%                   f1
              


[coord_1,~]=gll_arclengthmethod(N,f1,f1_prime,xa,ya,xb,yb);
[coord_2,~]=gll_arclengthmethod(N,f2,f2_prime,yc,xc,yb,xb);
[coord_3,~]=gll_arclengthmethod(N,f3,f3_prime,xd,yd,xc,yc);
[coord_4,~]=gll_arclengthmethod(N,f4,f4_prime,yd,xd,ya,xa);

[xi_gll,~]=gll(N);
[eta_gll,~]=gll(N);

X_xi=zeros(N-2,N-2);
Y_xi=zeros(N-2,N-2);

X_eta=zeros(N-2,N-2);
Y_eta=zeros(N-2,N-2);

X_eta_xi=zeros(N-2,N-2);
Y_eta_xi=zeros(N-2,N-2);

X_GH=zeros(N,N);
Y_GH=zeros(N,N);

X_GH(N,1:N)=coord_1(1:N,1);
X_GH(1:N,N)=coord_2(1:N,2);
X_GH(1,1:N)=coord_3(1:N,1);
X_GH(1:N,1)=coord_4(1:N,2);

Y_GH(N,1:N)=coord_1(1:N,2);
Y_GH(1:N,N)=coord_2(1:N,1);
Y_GH(1,1:N)=coord_3(1:N,2);
Y_GH(1:N,1)=coord_4(1:N,1);

phi_0= @(z) (1-z)/2;
phi_1= @(z) (1+z)/2;

for i=1:N-2
    for j=1:N-2
        X_xi(i,j) = phi_0(xi_gll(i+1))*X_GH(1,j+1)+phi_1(xi_gll(i+1))*X_GH(N,j+1);
        X_eta(i,j) = phi_0(eta_gll(j+1))*X_GH(i+1,1)+phi_1(eta_gll(j+1))*X_GH(i+1,N);
        X_eta_xi(i,j) = phi_0(xi_gll(i+1))*phi_0(eta_gll(j+1))*X_GH(1,1)...
                        +phi_1(xi_gll(i+1))*phi_0(eta_gll(j+1))*X_GH(N,1)...
                        +phi_1(xi_gll(i+1))*phi_1(eta_gll(j+1))*X_GH(N,N)...
                        +phi_0(xi_gll(i+1))*phi_1(eta_gll(j+1))*X_GH(1,N);
        
        Y_xi(i,j) = phi_0(xi_gll(i+1))*Y_GH(1,j+1)+phi_1(xi_gll(i+1))*Y_GH(N,j+1);
        Y_eta(i,j) = phi_0(eta_gll(j+1))*Y_GH(i+1,1)+phi_1(eta_gll(j+1))*Y_GH(i+1,N);
        Y_eta_xi(i,j) = phi_0(xi_gll(i+1))*phi_0(eta_gll(j+1))*Y_GH(1,1)...
                        +phi_1(xi_gll(i+1))*phi_0(eta_gll(j+1))*Y_GH(N,1)...
                        +phi_1(xi_gll(i+1))*phi_1(eta_gll(j+1))*Y_GH(N,N)...
                        +phi_0(xi_gll(i+1))*phi_1(eta_gll(j+1))*Y_GH(1,N);
                
        X_GH(i+1,j+1)=X_xi(i,j)+X_eta(i,j)-X_eta_xi(i,j);
        Y_GH(i+1,j+1)=Y_xi(i,j)+Y_eta(i,j)-Y_eta_xi(i,j);
    end
end
end



