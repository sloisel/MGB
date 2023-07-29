function [f0,f1,f2] = F(B,u,c)
%    [f0,f1,f2] = Barrier.F(u,c)
%
%    Input:
%    u    An (n×1) input.
%    c    (Optional). An (n×1) weight vector. The default value is c=0.
%
%    Outputs:
%    f0   The barrier value, calculated by:
%             f0 = dot(c,u) + dot(W,F0(D{1}*u,D{2}*u,...))
%    f1   The gradient of the barrier.
%    f2   The Hessian of the barrier.
%
%    If f1 and f2 are not requested, they will not be computed, potentially
%    saving some calculation time. f2 will be sparse if possible.
if(nargin<3), c=zeros(size(u)); end
d = length(B.D);
vv = cell(d,1);
n = size(u,1);
for pp=1:d
    vv{pp} = B.D{pp}*u;
end
m = size(vv{1},1);
f0 = dot(c,u) + dot(B.W,B.F0(vv{:}));
if(nargout<2), return; end
f1 = c;
for pp=1:d
    f1 = f1 + (B.D{pp})'*(B.W.*B.F1{pp}(vv{:}));
end
if(nargout<3), return; end
f2 = sparse(n,n);
for pp=1:d
    f2 = f2 + B.D{pp}'*spdiags(B.W.*B.F2{pp,pp}(vv{:}),0,m,m)*B.D{pp};
    for qq=1:pp-1
        Fpq = B.W.*B.F2{pp,qq}(vv{:});
        f2 = f2 + (B.D{pp}'*spdiags(Fpq,0,m,m)*B.D{qq}) ...
            +(B.D{qq}'*spdiags(Fpq,0,m,m)*B.D{pp});
    end
end
end

