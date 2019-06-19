classdef StereoR < StereoInterface
    %STEREO Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        R_centre
        R_sigma
        thres_stop_R
        
        delta
        angleMat
    end
    
    methods
        function obj = StereoR(p, q, R_centre, R_sigma, R_sigma_stop, delta)
            %STEREO Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@StereoInterface(p, q);
            obj.thres_stop_R = sqrt(3) * R_sigma_stop;
            
            % default angular error allowed in p and q measurements
            if isempty(delta)
                delta = 0;
            end
            obj.delta = delta;
            
            if isempty(R_centre) || isempty(R_sigma)
                % default R range to search
                obj.R_centre = [0, 0, 0];
                obj.R_sigma = pi;
            else
                obj.R_centre = R_centre;
                obj.R_sigma = R_sigma;
            end
        end
        
        
        function [obj] = setContext(obj, block)
            %SETCONTEXT Necessary computation before checking bounds
            % Given rotation at block centre, find angles(Rp, q)
            R  = block.aa2mat();
            Rp = (R * obj.p')';
            obj.angleMat = obj.angles(Rp, obj.q);
        end
        
        function [block] = updateLowerBound(obj, block, thres_stop) 
            %UPDATELOWERBOUND Set lower bound at stopping threshold
            assert(~isempty(obj.angleMat), 'Context was not set');
            block.edges_stop = obj.angleMat < thres_stop + 2*obj.delta;
            block.LB = obj.getMaxBipartiteMatching(block.edges_stop);
        end
        
        function [block] = updateUpperBound(obj, block, thres_stop)
            %SETUPPERBOUND Set upper bound at sqrt(3) * sigma threshold
            assert(~isempty(obj.angleMat), 'Context was not set');
            if block.thres > thres_stop
                edges = obj.angleMat < block.thres + 2*obj.delta;
                block.UB = obj.getMaxBipartiteMatching(edges);
            else
                block.UB = block.LB;
            end
        end
        
        function [obj, solutions] = findSolutions(obj)
            init_blk = cube(obj.R_centre, obj.R_sigma, min(obj.Np, obj.Nq));
            obj = obj.bnb(init_blk.subdivide(), obj.thres_stop_R, false); 
            solutions = obj.solutions;
        end
        
    end
end

