function [ xnew ] = perfect_shuffle( p , q, x )
% Splits x into p stacks of q elements and deals one from each stack q times.
% length(x) = p*q.

  counter = 1;
  for jj = 1: q
     for ii = 1: p
        ndx           = (ii - 1) * q + jj;
        xnew(counter) = x(ndx);
        counter       = counter + 1;
     end
  end

end


