classdef (Abstract) FEM
    %FEM Abstract base class for finite element implementation.
    %   This is the specification of what a FEM class should implement in
    %   order to be used with Algorithm MGB.
    %   Properties:
    %   n        The number of vertices in the mesh, counting
    %            multiplicities. For a piecewise quadratic mesh with m
    %            elements in 2d, n = 6m.
    %   h        The average mesh edge length.
    %   alpha    The order of accuracy of the mesh. For piecewise
    %            quadratic, alpha=2.
    %   d        The dimension of the mesh, e.g. d=2 for two dimensions.
    %
    %   Methods:
    %   continuous(M):
    %       Produce basis functions for the space of continuous functions
    %       of the highest order supported by the FEM implementation.
    %   linear(M):
    %       Produce basis functions for the space of discontinuous
    %       functions whose order is one less than the highest order
    %       supported by the FEM implementation. For piecewise quadratic
    %       elements, this would be the space of piecewise linear
    %       discontinuous functions.
    %   operators(M):
    %       Compute the differential operators associated with the mesh M.
    %   subdivide(M1,C):
    %       Subdivide a mesh by bisecting edges.
    %   vertices(M):
    %       Return an array V containing the vertices of the mesh.
    
    properties (Abstract)
        n, h, alpha, d
    end
    
    methods (Abstract)
        continuous(M)
        linear(M)
        operators(M)
        subdivide(M1,C)
        vertices(M)
    end
end

