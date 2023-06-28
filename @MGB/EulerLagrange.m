function [uh,fh,f,F] = EulerLagrange(M,L,u,randomize)
%EULERLAGRANGE Compute data for Euler-Lagrange problem
%   [uh,fh,f,F] = EulerLagrange(M,L,u,randomize)
%   Inputs:
%   M           A mesh.
%   L           A Lagrangian, given as a function handle.
%   u           A chosen solution, given as a function handle.
%   randomize   An optional vertex randomization. The default is 0. If this
%               number is given, vertices of the mesh will be perturbed by
%               randomize*randn(). You should never need this unless you've
%               accidentally constructed a highly singular problem
%               instance.
%   Outputs:
%   uh          The solution u evaluated at vertices of the mesh.
%   fh          The corresponding forcing at vertices of the mesh.
%   f           The forcing as a function handle.
%   F           The forcing as a symbolic math expression.
if(nargin<4)
    randomize = 0; %1e-10;
end
d = M.mesh.d;
P = sym('p',[d 1]);
PC = num2cell(P);
L1 = L(PC{:});
E = gradient(L1,P);
X = sym('x',[d 1]);
XC = num2cell(X);
U = u(XC{:});
GU = gradient(U,X);
E2 = subs(E,P,GU);
F = divergence(E2,X);
f = matlabFunction(F,'Vars',X);
V = M.mesh.vertices();
V = V+randn(size(V))*randomize;
V = mat2cell(V',ones(1,size(V,2)));
uh = u(V{:})';
fh = f(V{:})';
end

