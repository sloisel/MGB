function [M2,R] = subdivide(M1,C)
%SUBDIVIDE Subdivide a piecewise quadratic triangular mesh
%   [M2,R] = subdivide(M1,C)
%
%   Inputs:
%   M1     A mesh to subdivide.
%
%   Outputs:
%   M2     The subdivided mesh, obtained by bisecting each edge.
%   R      A matrix whose columns represent the basis functions of the mesh
%          M1, as they interpolate down to the mesh M2.
X = M1.X;
Y = M1.Y;
R0=sparse([         0         0         0    1.0000         0         0
         0    0.5000   -0.1250    0.5000    0.2500   -0.1250
         0    1.0000         0         0         0         0
    0.3750         0         0    0.7500         0   -0.1250
    0.3750    0.7500   -0.1250         0         0         0
    1.0000         0         0         0         0         0
         0    1.0000         0         0         0         0
   -0.1250    0.5000         0    0.2500    0.5000   -0.1250
         0         0         0         0    1.0000         0
   -0.1250    0.7500    0.3750         0         0         0
         0         0    0.3750         0    0.7500   -0.1250
         0         0    1.0000         0         0         0
         0         0         0         0    1.0000         0
   -0.1250    0.2500   -0.1250    0.5000    0.5000         0
         0         0         0    1.0000         0         0
         0         0   -0.1250         0    0.7500    0.3750
   -0.1250         0         0    0.7500         0    0.3750
         0         0         0         0         0    1.0000
         0         0         0         0    1.0000         0
   -0.1250    0.2500   -0.1250    0.5000    0.5000         0
         0         0         0    1.0000         0         0
   -0.1250    0.5000         0    0.2500    0.5000   -0.1250
         0    0.5000   -0.1250    0.5000    0.2500   -0.1250
         0    1.0000         0         0         0         0]);

V0 = [0   0     1
      0   0.5   0.5
      0   1     0
      0.5 0     0.5
      0.5 0.5   0
      1   0     0];
if(size(X,1)==3)
    U = V0*X;
    V = V0*Y;
    foo = repmat({sparse(fliplr(V0))},1,size(X,2));
    R = blkdiag(foo{:});
    M2 = FEM2d(U,V);
    return;
end

% B = {@(x,y) 2.*(0.5-x-y).*(1-x-y)
%     @(x,y) 4.*y.*(1-x-y)
%     @(x,y) 2.*y.*(y-0.5)
%     @(x,y) 4.*x.*(1-x-y)
%     @(x,y) 4.*x.*y
%     @(x,y) 2.*x.*(x-0.5)};
T0 = [1 2 4
      3 5 2
      6 4 5
      2 4 5];
U = cell(4,1);
V = cell(4,1);
for j=1:4
    Xj = X(T0(j,:),:);
    Yj = Y(T0(j,:),:);
    U{j} = V0*Xj;
    V{j} = V0*Yj;
end
U = reshape(vertcat(U{:}),6,4*size(X,2));
V = reshape(vertcat(V{:}),6,4*size(X,2));
foo = repmat({R0},1,size(X,2));
R = blkdiag(foo{:});
%R = zeros(24,6);
%for k=1:6
%    B{k}(X(k,1),Y(k,1))
%    R(:,k) = B{k}(U(:),V(:));
%end
M2 = FEM2d(U,V);
if(nargin>=2)
    R = reshape(R*reshape(C,6*size(X,2),1),6,size(U,2));
end
end

