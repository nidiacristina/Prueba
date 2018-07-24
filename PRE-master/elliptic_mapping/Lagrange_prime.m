function yint = Lagrange_prime(x,y,xx)

n = length(x);
if length(y)~=n, error('x and y must be same length'); end
s = 0;
for i = 1:n
  product = 1;
  denominator=1;
  sum=0;
  for k = 1:n
    if i ~= k
      denominator = denominator/(x(i)-x(k));
      for l = 1:n
          if k ~= l
              product = product*(xx-x(l));
          end 
      end
      sum = sum + product;
    end
  end
  s = s+y(i)*sum/denominator;
end
yint = s;