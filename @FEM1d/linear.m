function [L] = linear(M)
foo = repmat({sparse([1;1])},size(M.X,2),1);
L = blkdiag(foo{:});
end

