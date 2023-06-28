classdef FEM2d < FEM
    %FEM2d Piecewise quadratic finite elements
    %   This class does piecewise quadratic discontinuous finite elements
    %   in 2d, for use in Algorithm MGB. One can implement finite elements
    %   of various orders and dimensions by mimicking the interface
    %   presented here.
    %
    %   Properties:
    %   X     (6×m) array of X coordinates of vertices of elements. There
    %         are n triangular elements, each containing 6 vertices. The
    %         vertices are arranged as follows:
    %         3
    %         |\
    %         2 5
    %         |  \
    %         1-4-6
    %   Y     (6×m) array of X coordinates of vertices of elements.
    %   n     n = 6*n is the number of vertices in the mesh, counting
    %         duplicates.
    %   h     The average edge length. Edge length is measured between
    %         corner vertices (i.e. vertices 1,3,6).
    %   alpha The order alpha = 2 of the basis functions.
    %   d     d=2 is the dimension of the basis functions.
    
    properties
        X, Y, n, h, alpha, d
    end
    
    methods
        function M = FEM2d(X,Y)
            %FEM Construct an instance of this class
            %   M = FEM(X,Y) creates a mesh with the given vertices.
            M.X = X;
            M.Y = Y;
            if(size(X,1)==3)
                foo = M.subdivide();
                M.X = foo.X;
                M.Y = foo.Y;
            end
            M.n = numel(M.X);
            mynorm = @(v) mean(sqrt(v(:,1)^2+v(:,2)^2));
            M.h = (mynorm([M.X(3,1)-M.X(1,1),M.Y(3,1)-M.Y(1,1)]) + ...
                mynorm([M.X(3,1)-M.X(6,1),M.Y(3,1)-M.Y(6,1)]) + ...
                mynorm([M.X(6,1)-M.X(1,1),M.Y(6,1)-M.Y(1,1)]))/3;
            M.alpha = 2;
            M.d = 2;
        end
    end
end

