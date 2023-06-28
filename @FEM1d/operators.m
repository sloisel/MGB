function [D,W] = operators(M)
X = M.X;
Rx = sparse([-1 1;-1 1]);
W0 = [0.5;0.5];
U = X(2,:)-X(1,:);
detA = U;
Ainv = 1./detA;

Gx = cell(size(X,2),1);
W = cell(size(X,2),1);
for k=1:size(X,2)
    Gx{k} = (Ainv(k)*Rx);
    W{k} = W0*abs(detA(k));
end
Dx=blkdiag(Gx{:});
D = {Dx};
W = vertcat(W{:});
end

