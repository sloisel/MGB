function patch(M,C)
%PATCH A helper function that calls patch(...)
%   This is a thin wrapper for patch(...), to draw a given piecewise
%   quadratic mesh. See fem.subdivide for some examples.
X = M.X;
T0 = [     1     2     4
     3     5     2
     6     4     5
     2     4     5];
C = reshape(C,size(X,1),size(X,2));
plot(X,C,'-k');
end

