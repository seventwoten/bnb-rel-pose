classdef block
    %BLOCK Represents a 2D/3D region for use in branch-and-bound
    
    properties
        centre   % block centre
        sigma    % block half-length 
        parentUB % parent block upper bound
        
        UB       % block upper bound
        LB       % block lower bound
        
        thres       % threshold for comparison with stopping threshold
        capacityMat % capacity matrix for use in maxflow
        edges_stop  % edges matrix of possible matches at stopping threshold
    end
    
    methods
        function obj = block(c, s, p)
            %BLOCK Construct an instance of this class
            if nargin ~= 0
                obj.centre = c;
                obj.sigma = s;
                obj.parentUB = p;
            end
        end
        
        function [out] = getMaxBipartiteMatching(obj, edges)
            %GETMAXBIPARTITEMATCHING Return size of max bipartite matching
            Np = size(edges, 1);
            Nq = size(edges, 2);
            obj = obj.setCapacityMatrix(Np, Nq, edges);
            
            % Call maxflow to get max bipartite matching
            G = digraph(obj.capacityMat);
            out = maxflow(G, 1, Np+Nq+2);
        end
        
        function [obj] = setCapacityMatrix(obj, Np, Nq, edges)
            %SETCAPACITYMATRIX Set capacityMat attribute to template matrix
            C = obj.getTemplateCapacityMatrix(Np, Nq);
            C(2:Np+1, 1+Np+1 : 1+Np+Nq) = edges;
            obj.capacityMat = C;
        end
        
    end
    
    methods (Static)
        function [out] = getTemplateCapacityMatrix(Np, Nq)
            % Prepare shared template capacity matrix
            persistent C;
            if isempty(C)
                n_C = 1 + Np + Nq + 1;
                C = zeros(n_C, n_C);
                C(1, 2:Np+1) = 1; % set edges from s to p nodes
                C(1+Np+1 : 1+Np+Nq, end) = 1; % set edges from q nodes to t
            end
            out = C;
        end
      
    end
end

