classdef FEM1d < FEM
    %FEM1d Piecewise linear finite elements in 1d
    %  This is less well-tested than FEM2d, so see
    %  there for better documentation.
    
    properties
        X, n, h, alpha, d
    end
    
    methods
        function M = FEM1d(X)
            %FEM Construct an instance of this class
            %   M = FEM(X,Y) creates a mesh with the given vertices.
            M.X = X;
            M.n = numel(M.X);
            M.h = max(abs(X(2,:)-X(1,:)));
            M.alpha = 1;
            M.d = 1;
        end
    end
end

