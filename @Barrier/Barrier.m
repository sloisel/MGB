classdef Barrier
    % BARRIER An object to represent barrier functions.
    %    B = Barrier(F,D,W)
    %
    %    Inputs:
    %    F    A Matlab function, e.g. F = @(x) -log(x). This function
    %         should be usable in vector form, e.g. F = @(x) x.^2 is good,
    %         but F = @(x) x^2 is bad. This function will be called with
    %         symbols (see doc syms) and then symbolically differentiated.
    %         The function F gets stored as B.F0, and its derivatives are
    %         stored as B.F1 and B.F2.
    %    D    A cell array of matrices, e.g. D = {[3;2]}. Each D{k} matrix
    %         gives rise to a corresponding xk inputs to the F function.
    %    W    An (n×1) weight vector describing how to sum F values.
    %
    %    Output:
    %    B    an object of class Barrier.
    %
    %    Methods:
    %    Barrier.F
    %
    %    Properties:
    %    v    an array of automatically generated symbols used for input
    %         names for the function Barrier.f0 to compute its symbolic
    %         derivatives.
    %    D    A cell array of matrices (e.g. numerical differentiation).
    %    W    An (n×1) weight vector describing how to sum F values.
    %    F0   A copy of the function F passed to the constructor.
    %    F1   The gradient of F0.
    %    F2   The Hessian of F0.
    properties
        v
        D
        W
        F0
        F1
        F2
    end
    methods
        function B = Barrier(F,D,W)
            B.D = D;
            B.W = W;
            d = length(D);
            v = sym('x',[1,d]);
            B.v = v;
            v1 = cell(1,length(v));
            for k=1:d
                v1{k} = v(k);
            end
            B.F0 = F;
            FF = F(v1{:});
            F1S = cell(d,1);
            F2S = cell(d,d);
            B.F1 = cell(d,1);
            B.F2 = cell(d,d);
            WW = spdiags(W,0,length(W),length(W));
            for j=1:d
                F1S{j} = diff(FF,v(j));
                B.F1{j} = matlabFunction(F1S{j},'Vars',v);
                for k=1:j
                    F2S{j,k} = diff(F1S{j},v(k));
                    F2S{k,j} = F2S{j,k};
                    B.F2{j,k} = matlabFunction(F2S{j,k},'Vars',v);
                    B.F2{k,j} = B.F2{j,k};
                end
            end
            %            B = @fembarrier;
        end
    end
end
