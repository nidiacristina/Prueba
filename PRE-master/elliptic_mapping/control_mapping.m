function [s,t]=control_mapping(N, type_mapping)

% type_mapping='BOUNDARY_ORTHOGONALITY'

switch type_mapping
    case 'LAPLACE'
        
        f1= @(x) 0;
        f2= @(x) 0;
        f3= @(x) 0;
        f4= @(x) 0;

        f1_prime= @(x) 0;
        f2_prime= @(x) 0;
        f3_prime= @(x) 0;
        f4_prime= @(x) 0;

        xa=0;
        ya=0;
        xb=1;
        yb=0;
        xc=1;
        yc=1;
        xd=0;
        yd=1;

        [s,t]=gordon_hall_mapping(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
        f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd);
    
    case 'ARC_LENGHT'
        
        [x,~]=gll(N);
        x=(1.+x)/2;

        
%         x=linspace(0,1,N);
        
        s=zeros(N,N);
        t=zeros(N,N);
        
        for i=1:N
            s(N,i)=1;
            t(i,N)=1;
        end
        
        s(1:N,1)=x(1:N);
        s(1:N,N)=x(1:N);
        t(1,1:N)=x(1:N);
        t(N,1:N)=x(1:N);
        
        for i=1:N
            for j=1:N
                A=eye(2,2)-([s(i,N) 0; 0 t(N,j)]-[s(i,1) 0; 0 t(1,j)])*[0 1;1 0];
                b=[s(i,1);t(1,j)];
                x=A\b;
                s(i,j)=x(1);
                t(i,j)=x(2);
            end
        end
        
    case 'BOUNDARY_ORTHOGONALITY'
        [x,~]=gll(N);
        [~,D1,D2]=baryc(x);
        
        s=zeros(N,N);
        t=zeros(N,N);
        
        I=eye(N,N);
        D1_xi=kron(I,D1);
        D1_eta=kron(D1,I);
        D2_xi=kron(I,D2);
        D2_eta=kron(D2,I);
        
        L_s=D2_xi+D2_eta;
        L_t=D2_xi+D2_eta;
        
        b_neumann_t=[2:N-1 , (N^2-N+2):N^2-1];
        b_dirichlet_t=[1:N:N^2, N:N:N^2];
        
        b_dirichlet_s=[1:N , (N^2-N+1):N^2];
        b_neumann_s=[(N+1):N:(N^2-N), (2*N):N:(N^2-N)];
        
        L_s(sort(b_dirichlet_s),:)=zeros(2*N,N^2);
        L_s(sort(b_dirichlet_s),sort(b_dirichlet_s))=eye(2*N);
        
        L_t(sort(b_dirichlet_t),:)=zeros(2*N,N^2);
        L_t(sort(b_dirichlet_t),sort(b_dirichlet_t))=eye(2*N);
        
        rhs_s=zeros(N^2,1);
        rhs_t=zeros(N^2,1);
        rhs_s(b_dirichlet_s(N+1:2*N))=1;
        rhs_t(b_dirichlet_t(N+1:2*N))=1;
        
        L_s(sort(b_neumann_s),:)=D1_eta(sort(b_neumann_s),:);
        L_t(sort(b_neumann_t),:)=D1_xi(sort(b_neumann_t),:);
        
        s_vect=L_s\rhs_s;
        t_vect=L_t\rhs_t;
        
        s=vector_to_matrix_lex(s_vect,N,N);
        t=vector_to_matrix_lex(t_vect,N,N);
         
        H0= @(x) (1+2*x)*(1-x)^2;
        H1= @(x) (3-2*x)*x^2;
        
        for i=2:N-1
            for j=2:N-1
                
                f_s= @(y) s(i,1)*H0(t(1,j)*H0(y)+t(N,j)*H1(y))+...
                s(i,N)*H1(t(1,j)*H0(y)+t(N,j)*H1(y))-y;
            
                f_t= @(y) t(1,j)*H0(s(i,1)*H0(y)+s(i,N)*H1(y))+...
                t(N,j)*H1(s(i,1)*H0(y)+s(i,N)*H1(y))-y;
                
                s(i,j)=fzero(f_s,[-0.2,1.2]);
                t(i,j)=fzero(f_t,[-0.2,1.2]);
            
            end
        end
        
        
        
end
        
        