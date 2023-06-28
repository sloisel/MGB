function [x,SOL] = minimize(B,x,c,R,varargin)
%MINIMIZE Minimize the function dot(c,x) + B.F(x) in the span of R
%   [x,SOL] = minimize(B,x,c,R)
%
%   Inputs:
%   B     A Barrier object.
%   x     An (n×1) input.
%   c     An (n×1) weight vector.
%   R     An (n×m) matrix of basis vectors.
%
%   Outputs:
%   x     The converged minimizer of dot(c,x) + B.F(x).
%         The difference dx between the input value of x and the output
%         value of x, will be in the column span of R.
%   SOL   Details of the iteration. SOL.converged should be true if it
%         worked, and false if there was a failure (typically, exceeded the
%         maximum iteration count).
%
%   The algorithm is a damped Newton iteration with backtracking line
%   search.
%
%   Example:
%   B = Barrier(@(x) x.^2,{[2;3]},[1;1]);
%   [x,SOL] = B.minimize(3,1,1)
SOL = {};
SOL.maxit = 100;
SOL.alpha = 0.01;
SOL.beta = 0.25;
SOL.tol = 0.1/sqrt(size(x,1));
SOL.converged = false;
SOL.fk = [];
SOL.NewtonDecrement = [];
SOL.stepsize = [];
SOL.x = x;
SOL.xs = x;
SOL.saveiterations = false;
SOL.timeout = 300;
try
    t = toc;
catch
    tic;
end    
for k=1:2:length(varargin)
    switch(varargin{k})
        case 'maxit'
            SOL.maxit = varargin{k+1};
        case 'timeout'
            SOL.timeout = varargin{k+1};
        otherwise
            error('Unknown argument');
    end
end

for k=1:SOL.maxit
    if(toc>SOL.timeout), error('Timeout during Newton step'); end
    [f0,f1,f2] = B.F(SOL.x,c);
    if(~utils.isreal(f0)), error('f0 must be real!'); end
    if(any(~utils.isreal(f1))), error('f1 must be real!'); end
    SOL.fk(k) = f0;
    f2 = R'*f2*R;
    f2 = f2 + speye(size(f2,1))*norm(f2,inf)*1e-15;
    n = R*(f2\(R'*f1));
    if(any(~utils.isreal(n))), keyboard; error('Newton step must be real!'); end
    g = dot(f1,n);
    SOL.NewtonDecrement(k) = sqrt(g);
    s = 1;
    while(true)
        xa = SOL.x - s*n;
        if(all(xa==SOL.x) || s<1e-16), SOL.converged = true; break; end
        fa = B.F(xa,c);
        if(utils.isreal(fa) && fa<=f0-s*SOL.alpha*g)
            break;
        end
        s = s*SOL.beta;
    end
    SOL.stepsize(k) = s;
    SOL.x = xa;
    x = xa;
    if(SOL.saveiterations)
        SOL.xs(:,k+1) = SOL.x;
    end
    if(SOL.NewtonDecrement(k)<SOL.tol), SOL.converged = true; end
    if(SOL.converged), break; end
end
end

