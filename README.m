%%
%   Algorithm MGB is described in the paper:
%   @article{loisel2023algorithm,
%     title={Algorithm {MGB} to solve highly nonlinear elliptic {PDEs} in 
%            $\tilde{O}(n)$ {FLOPS}},
%     author={Loisel, S{\'e}bastien},
%     journal={arXiv preprint arXiv:2306.10183},
%     year={2023}
%     }
%
%   Included in this package:
%
%   Barrier    An object to implement barrier functions.
%   FEM        An abstract base class that describes the interface for FEM
%              modules.
%   FEM2d      A piecewise quadratic FEM implementation in 2d.
%   MGB        An object that implements Algorithm MGB.
%   utils      A small package of utility functions.
%
%   To use this package, proceed as follows.
%   1. Create a mesh using the FEM2d object.
%   2. Create a corresponding MGB object, this will create the multigrid
%      hierarchy.
%   3. Create a Barrier object.
%   4. Create initial value for uh, and a forcing vector fh.
%   5. Solve the Euler-Lagrange problem by calling MGB.Minimize
%   6. Use MGB.mesh.patch() to view the solution.
%
%   The rest of this file is a runnable example.

X = [0 0
     0 1
     1 1];
Y = [0 1
     1 1
     0 0];
mesh = FEM2d(X,Y);
M = MGB(mesh,5);
f = @(ux,uy,s) -log(s.^(2) - ux.^2 - uy.^2) - 2*log(s);
B = Barrier(f,M.D(end,:),M.W);
uh = double((M.mesh.X+0.1*M.mesh.Y<0.2)|(M.mesh.X-M.mesh.Y>0.8));
fh = double(zeros(size(M.mesh.X,1),size(M.mesh.X,2)));
[x,SOL] = M.minimize(uh,fh,B);
M.mesh.patch(x);
