function [M2,R] = subdivide(M1,C)
X = M1.X;
M = mean(X,1);
M2 = FEM1d(reshape([X(1,:);M;M;X(2,:)],2,size(X,2)*2));
R0= sparse([1   0
            0.5 0.5
            0.5 0.5
            0   1]);
R1 = repmat({R0},size(X,2),1);
R = blkdiag(R1{:});

if(nargin>=2)
    R = reshape(R*reshape(C,2*size(X,2),1),2,size(M2.X,2));
end
end

