classdef MGB
    %MGB This class implements Algorithm MGB.
    
    properties
        mesh
        D
        W
        Rs
        Cs
        CB
        Ls
        I
        derivatives
    end
    
    methods
        function M = MGB(mesh,L)
            Rs = {};
            Cs = {};
            CB = {};
            [Cs{1},CB{1}] = mesh.continuous();
            Ls = {mesh.linear()};
            for l=1:L-1
                [mesh,R] = mesh.subdivide();
                Rs{l} = R;
                [Cs{l+1},CB{l+1}] = mesh.continuous();
                Ls{l+1} = mesh.linear();
            end
            Rs{L} = speye(mesh.n);
            for l = L-1:-1:1
                Rs{l} = Rs{l+1}*Rs{l};
                Cs{l} = Rs{l}*Cs{l};
                CB{l} = Rs{l}*CB{l};
                Ls{l} = Rs{l}*Ls{l};
            end
            [D,M.W] = mesh.operators();
            Z = sparse(size(D{1},1),size(D{1},2));
            for k=1:length(D)
                M.D{k} = [D{k},Z];
            end
            M.D{end+1} = [Z,speye(size(Z,1))];
            for l=1:L
                M.I{l,1} = blkdiag(Cs{l},Ls{l});
            end
            M.mesh = mesh;
            M.Rs = Rs;
            M.Cs = Cs;
            M.CB = CB;
            M.Ls = Ls;
            M.derivatives = D;
        end
    end
end

