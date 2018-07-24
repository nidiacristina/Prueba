function L=arclength(f_prime,xa,xb)
%calculate the arclenght of function f between xa and xb

L=integral(@(x) (sqrt(1 + f_prime(x)^2)),xa,xb,'ArrayValued',true); 
