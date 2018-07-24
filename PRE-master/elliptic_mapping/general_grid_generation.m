function [X_grid,Y_grid]=general_grid_generation(N_1,M_1,N_2,M_2,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd,general_method,subdomain_method)

switch general_method
    case 'GH'
        n_iter_max_1=0;
        [X_gen,Y_gen]=elliptic_mapping(N_1,M_1,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
        f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd,n_iter_max_1);

    case 'elliptic'
        n_iter_max_1=1;
        [X_gen,Y_gen]=elliptic_mapping(N_1,M_1,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
        f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd,n_iter_max_1);

    case 'equidistant'
        n_iter_max_1=0;
        [X_gen,Y_gen]=elliptic_mapping(N_1,M_1,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
        f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd,n_iter_max_1);
    
        X_compu=linspace(-1,1,N_1);
        Y_compu=linspace(-1,1,M_1);
        
        [X_ini,~]=gll(N_1);
        [Y_ini,~]=gll(M_1);
        
%         X_ini=X_ini+1./2;
%         Y_ini=Y_ini+1./2;
        
        M_x=vandermonde_x(X_ini, X_compu,N_1);
        M_y=vandermonde_y(Y_ini, Y_compu,M_1);
        
        X_gen=M_x*X_gen'*M_y;
        Y_gen=M_x*Y_gen'*M_y;
        
end

figure
mesh(X_gen,Y_gen,zeros(24,8))
view(2);
axis('equal');

switch subdomain_method
    case 'GH'
        n_iter_max_2=0;
    case 'elliptic'
        n_iter_max_2=1;
end



M=M_1;
N=N_1;
%% First column 

f1_el= @(x) x*(Y_gen(2,2)-Y_gen(2,1))/(X_gen(2,2)-X_gen(2,1));
f2_el= @(x) x*(X_gen(1,2)-X_gen(2,2))/(Y_gen(1,2)-Y_gen(2,2));

f1_el_prime= @(x) (Y_gen(2,2)-Y_gen(2,1))/(X_gen(2,2)-X_gen(2,1));
f2_el_prime= @(x) (X_gen(1,2)-X_gen(2,2))/(Y_gen(1,2)-Y_gen(2,2));

[X_grid, Y_grid]=elliptic_mapping(N_2,M_2,f1_el,f1_el_prime,f2_el,f2_el_prime,f3,f3_prime,...
        f4,f4_prime,X_gen(2,1),Y_gen(2,1),X_gen(2,2),Y_gen(2,2),X_gen(1,2),Y_gen(1,2),...
        X_gen(1,1),Y_gen(1,1),n_iter_max_2);
    
for j=2:(M-2)
    
    f1_el= @(x) x*(Y_gen(j+1,2)-Y_gen(j+1,1))/(X_gen(j+1,2)-X_gen(j+1,1));
    f2_el= @(x) x*(X_gen(j,2)-X_gen(j+1,2))/(Y_gen(j,2)-Y_gen(j+1,2));
    f3_el= @(x) x*(Y_gen(j,1)-Y_gen(j,2))/(X_gen(j,1)-X_gen(j,2));

    f1_el_prime= @(x) (Y_gen(j+1,2)-Y_gen(j+1,1))/(X_gen(j+1,2)-X_gen(j+1,1));
    f2_el_prime= @(x) (X_gen(j,2)-X_gen(j+1,2))/(Y_gen(j,2)-Y_gen(j+1,2));
    f3_el_prime= @(x) (Y_gen(j,1)-Y_gen(j,2))/(X_gen(j,1)-X_gen(j,2));

    [X_el, Y_el]=elliptic_mapping(N_2,M_2,f1_el,f1_el_prime,f2_el,f2_el_prime,f3_el,f3_el_prime,...
        f4,f4_prime,X_gen(j+1,1),Y_gen(j+1,1),X_gen(j+1,2),Y_gen(j+1,2),X_gen(j,2),Y_gen(j,2),...
        X_gen(j,1),Y_gen(j,1),n_iter_max_2);
    
    X_grid=[X_grid; X_el];
    Y_grid=[Y_grid; Y_el];

end


f2_el= @(x) x*(X_gen(M-1,2)-X_gen(M,2))/(Y_gen(M-1,2)-Y_gen(M,2));
f3_el= @(x) x*(Y_gen(M-1,1)-Y_gen(M-1,2))/(X_gen(M-1,1)-X_gen(M-1,2));

f2_el_prime= @(x) (X_gen(M-1,2)-X_gen(M,2))/(Y_gen(M-1,2)-Y_gen(M,2));
f3_el_prime= @(x) (Y_gen(M-1,1)-Y_gen(M-1,2))/(X_gen(M-1,1)-X_gen(M-1,2));

[X_el, Y_el]=elliptic_mapping(N_2,M_2,f1,f1_prime,f2_el,f2_el_prime,f3_el,f3_el_prime,...
    f4,f4_prime,X_gen(M,1),Y_gen(M,1),X_gen(M,2),Y_gen(M,2),X_gen(M-1,2),Y_gen(M-1,2),...
    X_gen(M-1,1),Y_gen(M-1,1),n_iter_max_2);

X_grid=[X_grid; X_el];
Y_grid=[Y_grid; Y_el];





%%
for i=2:(N-2)

    f1_el= @(x) x*(Y_gen(2,i+1)-Y_gen(2,i))/(X_gen(2,i+1)-X_gen(2,i));
    f2_el= @(x) x*(X_gen(1,i+1)-X_gen(2,i+1))/(Y_gen(1,i+1)-Y_gen(2,i+1));
    f4_el= @(x) x*(X_gen(2,i)-X_gen(1,i))/(Y_gen(2,i)-Y_gen(1,i));

    f1_el_prime= @(x) (Y_gen(2,i+1)-Y_gen(2,i))/(X_gen(2,i+1)-X_gen(2,i));
    f2_el_prime= @(x) (X_gen(1,i+1)-X_gen(2,i+1))/(Y_gen(1,i+1)-Y_gen(2,i+1));
    f4_el_prime= @(x) (X_gen(2,i)-X_gen(1,i))/(Y_gen(2,i)-Y_gen(1,i));

    [X_col, Y_col]=elliptic_mapping(N_2,M_2,f1_el,f1_el_prime,f2_el,f2_el_prime,f3,f3_prime,...
            f4_el,f4_el_prime,X_gen(2,i),Y_gen(2,i),X_gen(2,i+1),Y_gen(2,i+1),X_gen(1,i+1),Y_gen(1,i+1),...
            X_gen(1,i),Y_gen(1,i),n_iter_max_2);

    for j=2:(M-2)

        f1_el= @(x) x*(Y_gen(j+1,i+1)-Y_gen(j+1,i))/(X_gen(j+1,i+1)-X_gen(j+1,i));
        f2_el= @(x) x*(X_gen(j,i+1)-X_gen(j+1,i+1))/(Y_gen(j,i+1)-Y_gen(j+1,i+1));
        f3_el= @(x) x*(Y_gen(j,i)-Y_gen(j,i+1))/(X_gen(j,i)-X_gen(j,i+1));
        f4_el= @(x) x*(X_gen(j+1,i)-X_gen(j,i))/(Y_gen(j+1,i)-Y_gen(j,i));

        f1_el_prime= @(x) (Y_gen(j+1,i+1)-Y_gen(j+1,i))/(X_gen(j+1,i+1)-X_gen(j+1,i));
        f2_el_prime= @(x) (X_gen(j,i+1)-X_gen(j+1,i+1))/(Y_gen(j,i+1)-Y_gen(j+1,i+1));
        f3_el_prime= @(x) (Y_gen(j,i)-Y_gen(j,i+1))/(X_gen(j,i)-X_gen(j,i+1));
        f4_el_prime= @(x) (X_gen(j+1,i)-X_gen(j,i))/(Y_gen(j+1,i)-Y_gen(j,i));

        [X_el, Y_el]=elliptic_mapping(N_2,M_2,f1_el,f1_el_prime,f2_el,f2_el_prime,f3_el,f3_el_prime,...
            f4_el,f4_el_prime,X_gen(j+1,i),Y_gen(j+1,i),X_gen(j+1,i+1),Y_gen(j+1,i+1),X_gen(j,i+1),Y_gen(j,i+1),...
            X_gen(j,i),Y_gen(j,i),n_iter_max_2);

        X_col=[X_col; X_el];
        Y_col=[Y_col; Y_el];

    end

    f1_el= @(x) f1(x+X_gen(M,i))-Y_gen(M,i);
    f2_el= @(x) x*(X_gen(M-1,i+1)-X_gen(M,i+1))/(Y_gen(M-1,i+1)-Y_gen(M,i+1));
    f3_el= @(x) x*(Y_gen(M-1,i)-Y_gen(M-1,i+1))/(X_gen(M-1,i)-X_gen(M-1,i+1));
    f4_el= @(x) x*(X_gen(M,i)-X_gen(M-1,i))/(Y_gen(M,i)-Y_gen(M-1,i));

    f1_el_prime= @(x) f1_prime(x+X_gen(M,i));
    f2_el_prime= @(x) (X_gen(M-1,i+1)-X_gen(M,i+1))/(Y_gen(M-1,i+1)-Y_gen(M,i+1));
    f3_el_prime= @(x) (Y_gen(M-1,i)-Y_gen(M-1,i+1))/(X_gen(M-1,i)-X_gen(M-1,i+1));
    f4_el_prime= @(x) (X_gen(M,i)-X_gen(M-1,i))/(Y_gen(M,i)-Y_gen(M-1,i));

    [X_el, Y_el]=elliptic_mapping(N_2,M_2,f1_el,f1_el_prime,f2_el,f2_el_prime,f3_el,f3_el_prime,...
        f4_el,f4_el_prime,X_gen(M,i),Y_gen(M,i),X_gen(M,i+1),Y_gen(M,i+1),X_gen(M-1,i+1),Y_gen(M-1,i+1),...
        X_gen(M-1,i),Y_gen(M-1,i),n_iter_max_2);

    X_col=[X_col; X_el];
    Y_col=[Y_col; Y_el];
    
    X_grid=[X_grid X_col];
    Y_grid=[Y_grid Y_col];
    
  
end
%%
%%Last column

f1_el= @(x) x*(Y_gen(2,N)-Y_gen(2,N-1))/(X_gen(2,N)-X_gen(2,N-1));
f4_el= @(x) x*(X_gen(2,N-1)-X_gen(1,N-1))/(Y_gen(2,N-1)-Y_gen(1,N-1));

f1_el_prime= @(x) (Y_gen(2,N)-Y_gen(2,N-1))/(X_gen(2,N)-X_gen(2,N-1));
f4_el_prime= @(x) (X_gen(2,N-1)-X_gen(1,N-1))/(Y_gen(2,N-1)-Y_gen(1,N-1));

[X_col, Y_col]=elliptic_mapping(N_2,M_2,f1_el,f1_el_prime,f2,f2_prime,f3,f3_prime,...
        f4_el,f4_el_prime,X_gen(2,N-1),Y_gen(2,N-1),X_gen(2,N),Y_gen(2,N),X_gen(1,N),Y_gen(1,N),...
        X_gen(1,N-1),Y_gen(1,N-1),n_iter_max_2);
    
for j=2:(M-2)
    
    f1_el= @(x) x*(Y_gen(j+1,N)-Y_gen(j+1,N-1))/(X_gen(j+1,N)-X_gen(j+1,N-1));
    f3_el= @(x) x*(Y_gen(j,N-1)-Y_gen(j,N))/(X_gen(j,N-1)-X_gen(j,N));
    f4_el= @(x) x*(X_gen(j+1,N-1)-X_gen(j,N-1))/(Y_gen(j+1,N-1)-Y_gen(j,N-1));

    f1_el_prime= @(x) (Y_gen(j+1,N)-Y_gen(j+1,N-1))/(X_gen(j+1,N)-X_gen(j+1,N-1));
    f3_el_prime= @(x) (Y_gen(j,N-1)-Y_gen(j,N))/(X_gen(j,N-1)-X_gen(j,N));
    f4_el_prime= @(x) (X_gen(j+1,N-1)-X_gen(j,N-1))/(Y_gen(j+1,N-1)-Y_gen(j,N-1));

    [X_el, Y_el]=elliptic_mapping(N_2,M_2,f1_el,f1_el_prime,f2,f2_prime,f3_el,f3_el_prime,...
        f4_el,f4_el_prime,X_gen(j+1,N-1),Y_gen(j+1,N-1),X_gen(j+1,N),Y_gen(j+1,N),X_gen(j,N),Y_gen(j,N),...
        X_gen(j,N-1),Y_gen(j,N-1),n_iter_max_2);
    
    X_col=[X_col; X_el];
    Y_col=[Y_col; Y_el];

end

f1_el= @(x) f1(x+X_gen(M,N-1))-Y_gen(M,N-1);
f3_el= @(x) x*(Y_gen(M-1,N-1)-Y_gen(M-1,N))/(X_gen(M-1,N-1)-X_gen(M-1,N));
f4_el= @(x) x*(X_gen(M,N-1)-X_gen(M-1,N-1))/(Y_gen(M,N-1)-Y_gen(M-1,N-1));

f1_el_prime= @(x) f1_prime(x+X_gen(M,N-1));
f3_el_prime= @(x) (Y_gen(M-1,N-1)-Y_gen(M-1,N))/(X_gen(M-1,N-1)-X_gen(M-1,N));
f4_el_prime= @(x) (X_gen(M,N-1)-X_gen(M-1,N-1))/(Y_gen(M,N-1)-Y_gen(M-1,N-1));


[X_el, Y_el]=elliptic_mapping(N_2,M_2,f1_el,f1_el_prime,f2,f2_prime,f3_el,f3_el_prime,...
    f4_el,f4_el_prime,X_gen(M,N-1),Y_gen(M,N-1),X_gen(M,N),Y_gen(M,N),X_gen(M-1,N),Y_gen(M-1,N),...
    X_gen(M-1,N-1),Y_gen(M-1,N-1),n_iter_max_2);

X_col=[X_col; X_el];
Y_col=[Y_col; Y_el];

X_grid=[X_grid X_col];
Y_grid=[Y_grid Y_col];



end

