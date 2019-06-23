classdef StereoT < StereoInterface
    %STEREOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        n1
        n2
        thres_stop_t
        t_long_lat
        t_half_len
        
        check_epipole  % option to turn on near-epipole check
        max_edges_p    % reject p points near epipole with edges > max_edges_p
        max_edges_q    % reject q points near epipole with edges > max_edges_q
        
        angleMat1
        angleMat2
    end
    
    methods
        function obj = StereoT(p, q, n1, n2, t_long_lat, t_half_len, t_half_len_stop, epipole_threshold)
            %STEREOT Construct an instance of this class
            %   t_long_lat and t_half_len are optional 
            %   (pass in [] to use default values)
            obj = obj@StereoInterface(p, q);
            obj.n1 = n1;
            obj.n2 = n2;
            obj.thres_stop_t = acos(cos(t_half_len_stop)^2);
            
            if isempty(t_long_lat) || isempty(t_half_len)
                % default t range to search
                obj.t_long_lat = [0, pi];
                obj.t_half_len = pi/2;
            else
                obj.t_long_lat = t_long_lat;
                obj.t_half_len = t_half_len;
            end
            
            if isempty(epipole_threshold)
                % default: reject those matching >80% of points
                epipole_threshold = 0.8;
            end
            
            if epipole_threshold == -1
                obj.check_epipole = false;
            else
                obj.check_epipole = true;
                obj.max_edges_p = epipole_threshold * obj.Nq;
                obj.max_edges_q = epipole_threshold * obj.Np;
            end
            
            
        end
        
        function [obj] = setContext(obj, block)
            obj.angleMat1 = obj.findAnglesToNormals(block.centre_xyz, obj.n1);
            obj.angleMat2 = obj.findAnglesToNormals(block.centre_xyz, obj.n2);
        end
        
        function [block] = updateLowerBound(obj, block, thres_stop) 
            %UPDATELOWERBOUND Update block lower bound at stopping threshold
            assert(~isempty(obj.angleMat1) & ~isempty(obj.angleMat2), 'Context was not set');
            positiveRange = pi/2 + thres_stop;
            block.edges_stop = ((obj.angleMat1 < positiveRange) & (obj.angleMat2 < positiveRange));
            
            % Near-epipole check: Remove rows/columns that match too many points
            if obj.check_epipole
                rows = sum(block.edges_stop,2) >= obj.max_edges_p;
                cols = sum(block.edges_stop,1) >= obj.max_edges_q;
                block.edges_stop(rows,:)= 0;
                block.edges_stop(:,cols)= 0;
            end
            
            block.LB = obj.getMaxBipartiteMatching(block.edges_stop);
        end
        
        function [block] = updateUpperBound(obj, block)
            %UPDATEUPPERBOUND Update block upper bound at threshold
            assert(~isempty(obj.angleMat1) & ~isempty(obj.angleMat2), 'Context was not set');
            positiveRange = pi/2 + block.thres;
            edges = ((obj.angleMat1 < positiveRange) & (obj.angleMat2 < positiveRange));
            block.UB = obj.getMaxBipartiteMatching(edges);
        end
        
        function [obj, solutions] = findSolutions(obj, early_stop)
            init_blk = sphericalPatch(obj.t_long_lat, obj.t_half_len);
            obj = obj.bnb(init_blk.subdivide(), obj.thres_stop_t, early_stop); 
            solutions = obj.solutions;
        end
        
        function [angleMat] = findAnglesToNormals(obj, t, n)
            %FINDANGLES Find angles between t-vector at patch centre, and normals in n
            [s1,s2,s3] = size(n);
            n = reshape(n, [], 3);
            angleMat = reshape(obj.angles(t, n), s1, s2);
        end
    end
end

