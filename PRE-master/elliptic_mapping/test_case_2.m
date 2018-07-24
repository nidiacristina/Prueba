clear all 
clc

N=25;


f1= @(x) sin(5/4*x)/2;
% f1= @(x) 0;
f2= @(x) sin(3/2*x)/2;
f3= @(x) 0;
f4= @(x) 0;

f1_prime= @(x) cos(5/4*x)*5/8;
% f1_prime= @(x) 0;
f2_prime= @(x) cos(3/2*x)*3/4;
f3_prime= @(x) 0;
f4_prime= @(x) 0;
X=zeros(24,1);
Y=zeros(24,1);


% for j=[1:2]
%     for i=[1:4]
%         fx= @(z) i*4*pi/5+f1(j*2*pi/3+f1(z))-z;
%         fy= @(z) j*2*pi/3+f1(i*4*pi/5+f1(z))-z;
%         x0=[i*4*pi/5-0.51 i*4*pi/5+0.51];
%         y0=[j*2*pi/3-0.51 j*2*pi/3+0.51];
%         X(j*6+i+1)=fzero(fx,x0);
%         Y(j*6+i+1)=j*2*pi/3+f1(X(j*6+i+1));
%         X(j*6+i+1)=i*4*pi/5+f1(Y(j*6+i+1));
%         
%     end
% end

for j=[0:3]
    for i=[0:5]
        X(j*6+i+1)=i*4*pi/5;
        Y(j*6+i+1)=j*2*pi/3;
    end
end
% 
% for j=[0:3]
%     for i=[0 5]
%         X(j*6+i+1)=i*4*pi/5;
%         Y(j*6+i+1)=j*2*pi/3;
%     end
% end



% figure
% hold on
% 
% fplot('4*pi/5+sin(x)/4',[0,2*pi])
% fplot('8*pi/5+sin(x)/4',[0,2*pi])
% fplot('12*pi/5+sin(x)/4',[0,2*pi])
% fplot('16*pi/5+sin(x)/4',[0,2*pi])
% 
% 
X_vect=vector_to_matrix_lex(X,4,6);
Y_vect=vector_to_matrix_lex(Y,4,6);
% 
plot(Y_vect,X_vect,'*')

% [X_GH,Y_GH]=gordon_hall_mapping_2(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
% f4,f4_prime,xa,ya,xb,yb,xc,yc,xd,yd);

[X_el_1,Y_el_1]=gordon_hall_mapping(N,f3,f3_prime,f2,f2_prime,f1,f1_prime,...
f3,f3_prime,X_vect(4,1),Y_vect(4,1),X_vect(4,2),Y_vect(4,2),X_vect(3,2),Y_vect(3,2),X_vect(3,1),Y_vect(3,1));
[X_el_2,Y_el_2]=elliptic_mapping(N,f1,f1_prime,f2,f2_prime,f1,f1_prime,...
f3,f3_prime,X(7),Y(7),X(8),Y(8),X(14),Y(14),X(13),Y(13));
[X_el_3,Y_el_3]=elliptic_mapping(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f3,f3_prime,X(13),Y(13),X(14),Y(14),X(20),Y(20),X(19),Y(19));
[X_el_4,Y_el_4]=elliptic_mapping(N,f3,f3_prime,f2,f2_prime,f1,f1_prime,...
f2,f2_prime,X(2),Y(2),X(3),Y(3),X(9),Y(9),X(8),Y(8));
[X_el_5,Y_el_5]=elliptic_mapping(N,f1,f1_prime,f2,f2_prime,f1,f1_prime,...
f2,f2_prime,X(8),Y(8),X(9),Y(9),X(15),Y(15),X(14),Y(14));
[X_el_6,Y_el_6]=elliptic_mapping(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f2,f2_prime,X(14),Y(14),X(15),Y(15),X(21),Y(21),X(20),Y(20));
[X_el_7,Y_el_7]=elliptic_mapping(N,f3,f3_prime,f2,f2_prime,f1,f1_prime,...
f2,f2_prime,X(3),Y(3),X(4),Y(4),X(10),Y(10),X(9),Y(9));
[X_el_8,Y_el_8]=elliptic_mapping(N,f1,f1_prime,f2,f2_prime,f1,f1_prime,...
f2,f2_prime,X(9),Y(9),X(10),Y(10),X(16),Y(16),X(15),Y(15));
[X_el_9,Y_el_9]=elliptic_mapping(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f2,f2_prime,X(15),Y(15),X(16),Y(16),X(22),Y(22),X(21),Y(21));
[X_el_10,Y_el_10]=elliptic_mapping(N,f3,f3_prime,f2,f2_prime,f1,f1_prime,...
f2,f2_prime,X(4),Y(4),X(5),Y(5),X(11),Y(11),X(10),Y(10));
[X_el_11,Y_el_11]=elliptic_mapping(N,f1,f1_prime,f2,f2_prime,f1,f1_prime,...
f2,f2_prime,X(10),Y(10),X(11),Y(11),X(17),Y(17),X(16),Y(16));
[X_el_12,Y_el_12]=elliptic_mapping(N,f1,f1_prime,f2,f2_prime,f3,f3_prime,...
f2,f2_prime,X(16),Y(16),X(17),Y(17),X(23),Y(23),X(22),Y(22));
[X_el_13,Y_el_13]=elliptic_mapping(N,f3,f3_prime,f3,f3_prime,f1,f1_prime,...
f2,f2_prime,X(5),Y(5),X(6),Y(6),X(12),Y(12),X(11),Y(11));
[X_el_14,Y_el_14]=elliptic_mapping(N,f1,f1_prime,f3,f3_prime,f1,f1_prime,...
f2,f2_prime,X(11),Y(11),X(12),Y(12),X(18),Y(18),X(17),Y(17));
[X_el_15,Y_el_15]=elliptic_mapping(N,f1,f1_prime,f3,f3_prime,f3,f3_prime,...
f2,f2_prime,X(17),Y(17),X(18),Y(18),X(24),Y(24),X(23),Y(23));

X_el=[X_el_3(1:end,1:end-1) X_el_6(1:end,1:end-1) X_el_9(1:end,1:end-1) X_el_12(1:end,1:end-1) X_el_15;...
    X_el_2(2:end,1:end-1) X_el_5(2:end,1:end-1) X_el_8(2:end,1:end-1) X_el_11(2:end,1:end-1) X_el_14(2:end,1:end);...
    X_el_1(2:end,1:end-1) X_el_4(2:end,1:end-1) X_el_7(2:end,1:end-1) X_el_10(2:end,1:end-1) X_el_13(2:end,1:end)];

Y_el=[Y_el_3(1:end,1:end-1) Y_el_6(1:end,1:end-1) Y_el_9(1:end,1:end-1) Y_el_12(1:end,1:end-1) Y_el_15;...
    Y_el_2(2:end,1:end-1) Y_el_5(2:end,1:end-1) Y_el_8(2:end,1:end-1) Y_el_11(2:end,1:end-1) Y_el_14(2:end,1:end);...
    Y_el_1(2:end,1:end-1) Y_el_4(2:end,1:end-1) Y_el_7(2:end,1:end-1) Y_el_10(2:end,1:end-1) Y_el_13(2:end,1:end)];

figure
mesh(X_el,Y_el,zeros(3*N-2,5*N-4))
% mesh(X_GH,Y_GH,zeros(N,N))
view(2);
axis('equal');
% 