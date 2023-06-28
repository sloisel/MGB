function [D,W] = operators(M)
%OPERATORS Differential operators for piecewise quadratic elements
%   [Dx,Dy,W] = operators(M)
%
%   Inputs:
%   M      A mesh.
%
%   Outputs:
%   D      D = {Dx,Dy}  
%          Matrices of the partial differential operators ∂/∂x and ∂/∂y.
%          If X,Y are (6×n) arrays then a coefficient vector u is expected
%          to be (6n×1), i.e. a column vector of dimension 6n. To
%          differentiate such a vector, use either Dx*u or Dy*u.
%          Note that if u has shape (6×n), it can be turned into a column
%          vector simply by doing reshape(u,6*n,1) or even u(:).
%   W      A (6n×1) vector so that dot(W,u) is the numerical integral of u.
%          The numerical quadrature is exact for piecewise quadratic
%          functions.

% Bx = {@(x,y) -3 + 4 * x + 4 * y
%       @(x,y) -4 * y
%       @(x,y) 0
%       @(x,y) 4 - 8 * x - 4 * y
%       @(x,y) 4 * y
%       @(x,y) 4 * x - 1};
% By = {@(x,y) -3 + 4 * x + 4 * y
%       @(x,y) 4 - 4 * x - 8 * y
%       @(x,y) 4 * y - 1
%       @(x,y) -4 * x
%       @(x,y) 4 * x
%       @(x,y) 0};
X = M.X;
Y = M.Y;
Rx = sparse([-3     0     0     4     0    -1
    -1    -2     0     2     2    -1
     1    -4     0     0     4    -1
    -1     0     0     0     0     1
     1    -2     0    -2     2     1
     1     0     0    -4     0     3]);
Ry = sparse([    -3     4    -1     0     0     0
    -1     0     1     0     0     0
     1    -4     3     0     0     0
    -1     2    -1    -2     2     0
     1    -2     1    -2     2     0
     1     0    -1    -4     4     0]);
W0 = [0;1;0;1;1;0]/6;
U = [X(6,:)-X(1,:);Y(6,:)-Y(1,:)];
V = [X(3,:)-X(1,:);Y(3,:)-Y(1,:)];
A = zeros(2,2,size(X,2));
A(:,1,:) = U;
A(:,2,:) = V;
detA = A(1,1,:).*A(2,2,:) - A(1,2,:).*A(2,1,:);
Ainv = zeros(2,2,size(X,2));
Ainv(1,1,:) = A(2,2,:)./detA;
Ainv(2,2,:) = A(1,1,:)./detA;
Ainv(1,2,:) = -A(1,2,:)./detA;
Ainv(2,1,:) = -A(2,1,:)./detA;

Gx = cell(size(X,2),1);
Gy = cell(size(X,2),1);
W = cell(size(X,2),1);
for k=1:size(X,2)
    Gx{k} = (Ainv(1,1,k)*Rx+Ainv(2,1,k)*Ry);
    Gy{k} = (Ainv(1,2,k)*Rx+Ainv(2,2,k)*Ry);
    W{k} = W0*abs(detA(k));
end
Dx=blkdiag(Gx{:});
Dy=blkdiag(Gy{:});
D = {Dx,Dy};
W = vertcat(W{:});
end

