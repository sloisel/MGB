function [L] = linear(M)
%LINEAR Calculate a basis for the piecewise linear discontinuous functions.
%   [L] = linear(M)
%   Output:
%   L    A matrix whose columns form a basis for the piecewise linear
%        discontinuous functions.
foo = repmat({sparse([1   0   0
               0.5 0.5 0
               0   1   0
               0.5 0   0.5
               0   0.5 0.5
               0   0   1])},size(M.X,2),1);
L = blkdiag(foo{:});
end

