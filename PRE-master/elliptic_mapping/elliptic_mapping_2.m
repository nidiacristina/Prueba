
function [X_el,Y_el]=elliptic_mapping_2(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd)

n_iter=0;
n_iter_max=1;

delta=0.05;

[x,~] = gll(N);
[~,D1,D2]=baryc(x);

[X0,Y0]=gordon_hall_mapping(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,f4,...
    f4_prime,xa,ya,xb,yb,xc,yc,xd,yd);
X1=zeros(N);
Y1=zeros(N);

X0_vect=matrix_to_vector_lex(X0);
Y0_vect=matrix_to_vector_lex(Y0);

X0_xi=vector_to_matrix_lex(kron(D1,eye(N))*X0_vect,N,N);
Y0_xi=vector_to_matrix_lex(kron(D1,eye(N))*Y0_vect,N,N);

X0_eta=vector_to_matrix_lex(kron(eye(N),D1)*X0_vect,N,N);
Y0_eta=vector_to_matrix_lex(kron(eye(N),D1)*Y0_vect,N,N);

X0_xi_xi=vector_to_matrix_lex(kron(D2,eye(N))*X0_vect,N,N);
Y0_xi_xi=vector_to_matrix_lex(kron(D2,eye(N))*Y0_vect,N,N);

X0_eta_eta=vector_to_matrix_lex(kron(eye(N),D2)*X0_vect,N,N);
Y0_eta_eta=vector_to_matrix_lex(kron(eye(N),D2)*Y0_vect,N,N);

X0_eta_xi=vector_to_matrix_lex(kron(D1,eye(N))*kron(eye(N),D1)*X0_vect,N,N);
Y0_eta_xi=vector_to_matrix_lex(kron(D1,eye(N))*kron(eye(N),D1)*Y0_vect,N,N);

S_0=zeros(N);
T_0=zeros(N);

for i=1:N
    for j=1:N
        g11=X0_xi(i,j)^2+Y0_xi(i,j)^2;
        g12=X0_eta(i,j)*X0_xi(i,j)+Y0_eta(i,j)*Y0_xi(i,j);
        g22=X0_eta(i,j)^2+Y0_eta(i,j)^2;
        R1=2*g12*X0_eta_xi(i,j)-g22*X0_xi_xi(i,j)-g11*X0_eta_eta(i,j);
        R2=2*g12*Y0_eta_xi(i,j)-g22*Y0_xi_xi(i,j)-g11*Y0_eta_eta(i,j);
        a=[g22*X0_xi(i,j) g11*X0_eta(i,j);g22*Y0_xi(i,j) g11*Y0_eta(i,j)]\[R1;R2];
        S_0(i,j)=a(1)*g22;
        T_0(i,j)=a(2)*g11;
    end
end

S_0_smooth=S_0;
T_0_smooth=T_0;


for i=2:N-1
    for j=2:N-1
%         S_0_smooth(i,j)=(S_0(i,j+1)+S_0(i,j-1))/2;
%         T_0_smooth(i,j)=(T_0(i+1,j)+T_0(i-1,j))/2;
        S_0(i,j)=(S_0(i,j+1)+S_0(i,j-1))/2;
        T_0(i,j)=(T_0(i+1,j)+T_0(i-1,j))/2;
    end
end

P=zeros(N);
Q=zeros(N);
R=zeros(N);
S_1=zeros(N);
T_1=zeros(N);
S=zeros(N);
T=zeros(N);

b=sort([1:N , (N+1):N:(N^2-N), (2*N):N:(N^2-N), (N^2-N+1):N^2]);
% figure
% hold on
% mesh(X0,Y0,zeros(N,N))
% view(2);
% axis('equal');
% 
% pause

while n_iter < n_iter_max
        
        
        X0_vect=matrix_to_vector_lex(X0);
        Y0_vect=matrix_to_vector_lex(Y0);

        X0_xi=vector_to_matrix_lex(kron(D1,eye(N))*X0_vect,N,N);
        Y0_xi=vector_to_matrix_lex(kron(D1,eye(N))*Y0_vect,N,N);

        X0_eta=vector_to_matrix_lex(kron(eye(N),D1)*X0_vect,N,N);
        Y0_eta=vector_to_matrix_lex(kron(eye(N),D1)*Y0_vect,N,N);
       
        X0_xi_xi=vector_to_matrix_lex(kron(D2,eye(N))*X0_vect,N,N);
        Y0_xi_xi=vector_to_matrix_lex(kron(D2,eye(N))*Y0_vect,N,N);

        X0_eta_eta=vector_to_matrix_lex(kron(eye(N),D2)*X0_vect,N,N);
        Y0_eta_eta=vector_to_matrix_lex(kron(eye(N),D2)*Y0_vect,N,N);

        for i=1:N
            for j=1:N
                P(i,j)=X0_eta(i,j)^2+Y0_eta(i,j)^2;
                Q(i,j)=X0_eta(i,j)*X0_xi(i,j)+Y0_eta(i,j)*Y0_xi(i,j);
                R(i,j)=X0_xi(i,j)^2+Y0_xi(i,j)^2;
                S_1(i,j)=-(X0_xi(i,j)*X0_xi_xi(i,j)+Y0_xi(i,j)*Y0_xi_xi(i,j))*P(i,j)/R(i,j)...
                    -(X0_xi(i,j)*X0_eta_eta(i,j)+Y0_xi(i,j)*Y0_eta_eta(i,j));
                T_1(i,j)=-(X0_eta(i,j)*X0_eta_eta(i,j)+Y0_eta(i,j)*Y0_eta_eta(i,j))/P(i,j)*R(i,j)...
                    -(X0_eta(i,j)*X0_xi_xi(i,j)+Y0_eta(i,j)*Y0_xi_xi(i,j));
                c=exp(-i*j*(N-i)*(N-j)/(delta*N^4));
                S(i,j)=c*S_0(i,j)+(1-c)*S_1(i,j);
                T(i,j)=c*T_0(i,j)+(1-c)*T_1(i,j);
            end
        end
        
        P_diag=diag(matrix_to_vector_lex(P));
        Q_diag=diag(matrix_to_vector_lex(Q));
        R_diag=diag(matrix_to_vector_lex(R));
        S_diag=diag(matrix_to_vector_lex(S));
        T_diag=diag(matrix_to_vector_lex(T));
        
        A=P_diag*kron(eye(N),D2)-2*Q_diag*kron(eye(N),D1)*kron(D1,eye(N))+...
            R_diag*kron(D2,eye(N))+ S_diag*kron(eye(N),D1)+T_diag*kron(D1,eye(N));
        
        A(b,:) = zeros(4*(N-1),N^2);
        A(b,b) = eye(4*(N-1));
        
        rhs_x=zeros(N^2,1);
        rhs_y=zeros(N^2,1);
        rhs_x(b)=X0_vect(b);
        rhs_y(b)=Y0_vect(b);
        
        X1=X0;
        Y1=Y0;
        
        X0=vector_to_matrix_lex(A\rhs_x,N,N);
        Y0=vector_to_matrix_lex(A\rhs_y,N,N);
        
        n_iter=n_iter+1;
        
%         if max(abs(X0-X1),abs(Y0-Y1))>10^(-1)
%             break
%         end



if n_iter==20
    disp('max iteration')
end;

X_el=X0;
Y_el=Y0;

% mesh(X_el,Y_el,zeros(N,N))
% 
% % mesh(X_GH,Y_GH,zeros(N,N))
% view(2);
% axis('equal');
% 
% pause
end
