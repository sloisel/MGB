function [CI,CB] = continuous(M)
%CONTINUOUS Basis for space of continuous functions
%   [CI,CB] = continuous(M)
%
%   Output:
%   CI     An (6n×k) matrix whose columns form a basis for the continuous
%          piecewise quadratic functions with homogeneous Dirichlet
%          conditions. This is determined by looking for vertices 
%          [X(i,j),Y(i,j)] that are repeated. Note that if
%          roundoff is present in the values of X and Y, then C might
%          break. One way to ensure that no roundoff is present in X and Y
%          is to start with integer coordinates before one applies
%          subdivision.
%   CB     Basis vectors for boundary values.
X = M.X;
Y = M.Y;
XY = [reshape(X,6*size(X,2),1) reshape(Y,6*size(Y,2),1)];
[C,iA,iC] = unique(XY,'rows');
R = sparse(1:size(XY,1),iC,1,size(XY,1),size(C,1));
T = reshape(iC,6,size(X,2));
E = [T([1,2],:),T([2,3],:),T([3,5],:),T([5 6],:),T([6 4],:),T([4 1],:)]';
E = sort(E,2);
[F,iA,iC] = unique(E,'rows');
F = accumarray(iC,1);
G = [E F(iC)];
idx = (G(:,3)==1);
B = G(idx,1:2); 
B = unique(B(:));
I = setdiff(1:size(R,2),B);
CI = R(:,I);
CB = R(:,B);
end

