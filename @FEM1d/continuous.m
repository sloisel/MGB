function [CI,CB] = continuous(M)
X = M.X;
foo = cell(size(X,2)+1,1);
I = sparse(1,1,1,1,1);
foo{1} = I;
foo{end} = I;
for k=2:size(X,2)
    foo{k} = sparse([1;1]);
end
R = blkdiag(foo{:});
CI = R(:,2:end-1);
CB = R(:,[1,size(X,2)+1]);
end

