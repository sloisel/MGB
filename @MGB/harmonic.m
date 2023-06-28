function [u] = harmonic(M,g)
%HARMONIC Compute the harmonic continuation of boundary data g.
%   [u] = harmonic(M,g)
%   Inputs:
%   M    A mesh.
%   g    A finite element function on M.
%   Outputs:
%   u    A finite element function on M. u and g will have the same
%        Dirichlet data, but u will be discretely harmonic.
n = size(M.derivatives{1},1);
L = sparse(n,n);
for k=1:length(M.derivatives)
    L = L + M.derivatives{k}'*spdiags(M.W,0,n,n)*M.derivatives{k};
end
CB = M.CB{end};
AII = M.Cs{end}'*L*M.Cs{end};
AIB = M.Cs{end}'*L*CB;
m = size(CB,2);
MM = spdiags(1./sum(CB,1)',0,m,m);
h = MM*(CB'*g);
%assert(norm(CB'*(CB*h)-CB'*g,inf)/norm(g,inf)<1e-10);
u = CB*h-(M.Cs{end}*(AII\(AIB*h)));
end

