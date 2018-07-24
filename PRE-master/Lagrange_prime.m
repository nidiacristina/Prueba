function yint = Lagrange_prime(x,y,xx)

n = length(x);
if length(y)~=n, error('x and y must be same length'); end
s = 0;
t = 0;
for i = 1:n
    for k=1:n
       if i ~= k
           product = y(i);
           for j = 1:n
             if i ~= j
               product = product*(xx-x(j))/(x(i)-x(j));
             end
           end
       end
       t=t+product/(xx-x(k));
    end
  s = s+t;
  t=0;
end
yint = s;