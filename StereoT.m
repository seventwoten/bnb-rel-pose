classdef StereoT < StereoInterface
    %STEREOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n1
        n2
        thres_stop_t
        t_long_lat
        t_half_len
        
        angleMat1
        angleMat2
    end
    
    methods
        function obj = StereoT(p, q, n1, n2, t_long_lat, t_half_len, thres_stop_t)
            %STEREOT Construct an instance of this class
            %   t_long_lat and t_half_len are optional 
            %   (pass in [] to use default values)
            obj = obj@StereoInterface(p, q);
            obj.n1 = n1;
            obj.n2 = n2;
            obj.thres_stop_t = thres_stop_t;
            
            if isempty(t_long_lat) || isempty(t_half_len)
                % default t range to search
                obj.t_long_lat = [0, pi];
                obj.t_half_len = pi/2;
            else
                obj.t_long_lat = t_long_lat;
                obj.t_half_len = t_half_len;
            end
        end
        
        function [obj] = setContext(obj, block)
            obj.angleMat1 = obj.angles(block.centre_xyz, obj.n1);
            obj.angleMat2 = obj.angles(block.centre_xyz, obj.n2);
        end
        
        function [block] = updateLowerBound(obj, block, thres_stop) 
            %UPDATELOWERBOUND Update block lower bound at stopping threshold
            assert(~isempty(obj.angleMat1) & ~isempty(obj.angleMat2), 'Context was not set');
            positiveRange = pi/2 + thres_stop;
            block.edges_stop = ((obj.angleMat1 < positiveRange) & (obj.angleMat2 < positiveRange));
            block.LB = obj.getMaxBipartiteMatching(block.edges_stop);
        end
        
        function [block] = updateUpperBound(obj, block, thres_stop)
            %UPDATEUPPERBOUND Update block upper bound at threshold
            assert(~isempty(obj.angleMat1) & ~isempty(obj.angleMat2), 'Context was not set');
            if block.thres > thres_stop
                positiveRange = pi/2 + block.thres;
                edges = ((obj.angleMat1 < positiveRange) & (obj.angleMat2 < positiveRange));
                block.UB = obj.getMaxBipartiteMatching(edges);
            else
                block.UB = block.LB;
            end
            
        end
        
        function [obj, solutions] = findSolutions(obj, early_stop)
            init_blk = sphericalPatch(obj.t_long_lat, obj.t_half_len, []);
            obj = obj.bnb(init_blk.subdivide(), obj.thres_stop_t, early_stop); 
            solutions = obj.solutions;
        end
        
    end
    
    methods (Static)
        function [angleMat] = angles(t, n)
            %ANGLES Find angles between t-vector at patch centre, and normals in n
            [s1,s2,s3] = size(n);
            n = reshape(n, [], 3);
            angleMat = reshape(angles(t, n), s1, s2);
        end
    end
end

