function [u,SOL] = minimize(M,uh,fh,B,varargin)
%MINIMIZE Solves the Euler-Lagrange problem with Algorithm MGB
%   Inputs:
%   M    An MGB object.
%   uh   The initial value for uh.
%   fh   The forcing.
%   B    a Barrier object.
%   Outputs:
%   u    The solution of the Euler-Lagrange problem.
%   SOL  A detailed solution object.
uh = uh(:);
fh = fh(:);
s0 = 1;
uh = M.harmonic(uh);
L = length(M.I);
while(true)
    s1 = s0*ones(size(M.Ls{end},1),1);
    x = [uh;s1];
    if(utils.isreal(B.F(x)))
        break;
    end
    s0 = abs(s0)*2+1;
end
c = [M.W.*fh
     M.W];
SOL.maxit = 100000;
SOL.tol = max(0.1*M.mesh.h^(M.mesh.alpha*2),1e-8);
SOL.c = c;
SOL.B = B;
SOL.rho = 2;
SOL.converged = false;
SOL.x = x;
SOL.xs = x;
SOL.M = M;
SOL.uh = uh;
SOL.fh = fh;
SOL.saveiterations = false;
SOL.varargin = varargin;
try
    SOL.t0 = toc;
catch
    tic;
end    
SOL.t0 = toc;
SOL.timeout = SOL.t0+300;
tinit = 0.1/size(uh,1);
SOL.ts = tinit;
SOL.ls = 1;
SOL.naive_maxit = 5;
SOL.mgb_maxit = SOL.maxit;
SOL.schedule = inf;
SOL.inner = cell(L+1,1);
SOL.its = zeros(L+1,1);
for k=1:2:length(varargin)
    switch(varargin{k})
        case 'schedule'
            SOL.schedule = varargin{k+1};
        case 'naive_maxit'
            SOL.naive_maxit = varargin{k+1};
        case 'mgb_maxit'
            SOL.mgb_maxit = varargin{k+1};
        otherwise
            error('Unknown argument');
    end
end
    function l = schedule(t)
        if(~isa(SOL.schedule,'function_handle'))
            l = SOL.schedule;
        else
            l = SOL.schedule(t);
        end
        l = round(min(L,max(1,l)));
    end
for k=2:SOL.maxit
    SOL.ts(k) = SOL.ts(k-1)*SOL.rho(k-1);
    SOL.ls(k) = schedule(SOL.ts(k));
    xprev = x;
    [x,SOL.inner{1,k}] = B.minimize(x,SOL.ts(k)*c,M.I{SOL.ls(k-1)},'maxit',SOL.naive_maxit,'timeout',SOL.timeout);
    SOL.its(1,k) = length(SOL.inner{1,k}.fk);
    if(~SOL.inner{1,k}.converged)
        if(SOL.mgb_maxit==0), error('Naïve step did not converge'); end
        x = xprev;
        for l = 1:SOL.ls(k-1)
            [x,SOL.inner{l+1,k}] = B.minimize(x,SOL.ts(k)*c,M.I{l},'maxit',SOL.mgb_maxit,'timeout',SOL.timeout);
            SOL.its(l+1,k) = length(SOL.inner{l+1,k}.fk);
            assert(SOL.inner{l+1,k}.converged);
        end
    end
    for l = SOL.ls(k-1)+1:SOL.ls(k)
        [x,SOL.inner{l+1,k}] = B.minimize(x,SOL.ts(k)*c,M.I{l},'maxit',SOL.maxit,'timeout',SOL.timeout);
        SOL.its(l+1,k) = length(SOL.inner{l+1,k}.fk);
        assert(SOL.inner{l+1,k}.converged);
    end
    m = max(SOL.its(:,k));
    if(m<=2)
        SOL.rho(k) = min(SOL.rho(k-1)^2,2);
    elseif(m>=6)
        SOL.rho(k) = sqrt(SOL.rho(k-1));
        if(SOL.rho(k)==1), error('t step size reduced to 1'); end
    else
        SOL.rho(k) = SOL.rho(k-1);
    end

    elapsed = toc-SOL.t0;
%    eta = (log(1/SOL.tol)-log(tinit))/(log(SOL.ts(k))-log(tinit))*elapsed;
%    remaining = eta-elapsed;
    disp(sprintf('Iteration %d, t=%f. Elapsed=%.1fs',k,SOL.ts(k),elapsed));

    if(SOL.ts(k)*SOL.tol>1), break; end
end
u = x(1:size(M.D{1},1),1);
end

