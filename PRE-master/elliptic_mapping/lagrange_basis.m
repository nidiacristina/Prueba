function y=lagrange_basis(X,m,xx)

n=length(X);
% s = 0;
% t = 0;
% for i = 1:n
%   product = 1;
%   for j = 1:n
%     if i ~= j
%       product = product*(X(i)-X(j));
%     end
%   end
%   
%   if i==m
%       t=1/product*(xx-X(i));
%   end
%   
%   s = s+1/product*(xx-X(i));
%   
% end
% y = t/s;
product=1;
for j = 1:n
    if m ~= j
      product = product*(xx-X(j))/(X(m)-X(j));
    end
end
y = product;
