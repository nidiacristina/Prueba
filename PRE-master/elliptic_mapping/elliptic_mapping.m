function [X_el,Y_el]=elliptic_mapping(N,M,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd,n_iter_max)

%Calculate the position of the points of the mapping with Gordon Hall
%algorithm
%                   f3
%       xd  x - - - - - - - x xc
%           |               |
%           |               |
%        f4 |               | f2 M
%           |               |
%           |               |
%           |               |
%       xa  x - - - - - - - x xb
%                   f1 N
%
%   n_iter_max=0 --> Gordon Hall
%   n_iter_max=1 --> Laplace

n_iter=0;


[x,~] = gll(N);
[~,D1_eta,D2_eta]=baryc(x);
[y,~] = gll(M);
[~,D1_xi,D2_xi]=baryc(y);

[X0,Y0]=gordon_hall_mapping(N,M,f1,f1_prime,f2,f2_prime,f3,f3_prime,f4,...
    f4_prime,xa,ya,xb,yb,xc,yc,xd,yd);

X1=zeros(M,N);
Y1=zeros(M,N);


P=zeros(M,N);
Q=zeros(M,N);
R=zeros(M,N);

b=sort([1:N , (N+1):N:(N*M-N), (2*N):N:(N*M-N), (N*M-N+1):N*M]);

while n_iter < n_iter_max
        
        
        X0_vect=matrix_to_vector_lex(X0);
        Y0_vect=matrix_to_vector_lex(Y0);

        X0_xi=vector_to_matrix_lex(kron(D1_xi,eye(N))*X0_vect,M,N);
        Y0_xi=vector_to_matrix_lex(kron(D1_xi,eye(N))*Y0_vect,M,N);

        X0_eta=vector_to_matrix_lex(kron(eye(M),D1_eta)*X0_vect,M,N);
        Y0_eta=vector_to_matrix_lex(kron(eye(M),D1_eta)*Y0_vect,M,N);
       
        for i=1:M
            for j=1:N
                P(i,j)=X0_eta(i,j)^2+Y0_eta(i,j)^2;
                Q(i,j)=X0_eta(i,j)*X0_xi(i,j)+Y0_eta(i,j)*Y0_xi(i,j);
                R(i,j)=X0_xi(i,j)^2+Y0_xi(i,j)^2;
            end
        end
        
        P_diag=diag(matrix_to_vector_lex(P));
        Q_diag=diag(matrix_to_vector_lex(Q));
        R_diag=diag(matrix_to_vector_lex(R));
        
        A=P_diag*kron(eye(M),D2_eta)-2*Q_diag*kron(eye(M),D1_eta)*kron(D1_xi,eye(N))+...
            R_diag*kron(D2_xi,eye(N));
        
        A(b,:) = zeros(2*(N-1)+2*(M-1),N*M);
        A(b,b) = eye(2*(N-1)+2*(M-1));
        
        rhs_x=zeros(N*M,1);
        rhs_y=zeros(N*M,1);
        rhs_x(b)=X0_vect(b);
        rhs_y(b)=Y0_vect(b);
        
        X1=X0;
        Y1=Y0;
        
        X0=vector_to_matrix_lex(A\rhs_x,M,N);
        Y0=vector_to_matrix_lex(A\rhs_y,M,N);
        
        n_iter=n_iter+1;
        
end
X_el=X0;
Y_el=Y0;


end





    
