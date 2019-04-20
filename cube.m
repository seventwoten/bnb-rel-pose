classdef cube < block
    %CUBE 3D block
    properties
        angleMat
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
            obj.edges_stop = obj.angleMat < thres_stop;
            obj.LB = obj.getMaxBipartiteMatching(obj.edges_stop);
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
        
    end
    
end

