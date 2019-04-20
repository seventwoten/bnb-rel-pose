classdef cube < block
    %CUBE 3D block
    properties
        thres
        angleMat
        capacityMat
        edges_tight
    end
    
    methods
        function obj = cube(c, s, p)
            %CUBE Construct an instance of this class
            if nargin == 0
                c = [];
                s = [];
                p = [];
            end
            obj = obj@block(c, s, p);
            obj.thres = sqrt(3) * obj.sigma;
        end
    
        function [subblocks] = subdivide(obj)
            %SUBDIVIDE Subdivide block and return 8 subblocks
            shifts = [-1, -1, -1 ;
                      -1, -1,  1 ;
                      -1,  1, -1 ;
                      -1,  1,  1 ;
                       1, -1, -1 ;
                       1, -1,  1 ;
                       1,  1, -1 ;
                       1,  1,  1  ];
            
            sigma_new  = obj.sigma/2;
            centre_new = obj.centre + sigma_new * shifts;
            
            for c = size(centre_new, 1):-1:1
                subblocks(c) = cube(centre_new(c,:), sigma_new, obj.UB);
            end
            
        end
        
        function [obj] = setContext(obj, p, q)
            %SETCONTEXT Necessary computation before checking bounds
            R  = aa2mat(obj.centre);
            Rp = (R * p')';
            obj.angleMat = angles(Rp, q);
        end
        
        function [obj] = setLowerBound(obj, thres_stop)
            %SETLOWERBOUND Set lower bound at stopping threshold
            assert(~isempty(obj.angleMat), 'Context was not set');
            obj.edges_tight = obj.angleMat < thres_stop;
            obj.LB = obj.getMaxBipartiteMatching(obj.edges_tight);
        end
        
        function [obj] = setUpperBound(obj, thres_stop)
            %SETUPPERBOUND Set upper bound at sqrt(3) * sigma threshold
            assert(~isempty(obj.angleMat), 'Context was not set');
            if obj.thres > thres_stop
                edges = obj.angleMat < obj.thres;
                obj.UB = obj.getMaxBipartiteMatching(edges);
            else
                obj.UB = obj.LB;
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

