classdef StereoRT < StereoInterface
    %STEREORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        R_centre
        R_sigma
        thres_stop_R
        
        t_long_lat
        t_half_len
        t_half_len_stop
        
        delta
        epipole_threshold
        
        n1_LB
        n2_LB
        n1_UB
        n2_UB
        angleMat
    end
    
    methods
        function obj = StereoRT(p, q, R_centre, R_sigma, R_sigma_stop, t_long_lat, t_half_len, t_half_len_stop, delta, epipole_threshold)
            %STEREORT Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@StereoInterface(p, q);
            obj.thres_stop_R = sqrt(3) * R_sigma_stop;
            obj.t_half_len_stop = t_half_len_stop;
            obj.t_long_lat = t_long_lat;
            obj.t_half_len = t_half_len;
            obj.delta = delta;
            obj.epipole_threshold = epipole_threshold;
            
            if isempty(R_centre) || isempty(R_sigma)
                % default t range to search
                obj.R_centre = [0, 0, 0];
                obj.R_sigma = pi;
            else
                obj.R_centre = R_centre;
                obj.R_sigma = R_sigma;
            end
        end
        
        function [obj] = setContext(obj, block)
            %SETCONTEXT Compute (n1, n2) from (R, p, q, angles, thres, delta)
            [obj.n1_LB, obj.n2_LB] = getWedges(obj, block, obj.thres_stop_R);
            
            % Skip upper bound computation if below stopping threshold
            if block.thres > obj.thres_stop_R 
                [obj.n1_UB, obj.n2_UB] = getWedges(obj, block, block.thres);
            end
        end
        
        function [normals1, normals2] = getWedges(obj, block, threshold_R)
            R  = block.aa2mat();
            Rp = (R * obj.p')';
            obj.angleMat = obj.angles(Rp, obj.q);
            
            sin_e_p = sin(obj.delta + threshold_R);
            sin_e_q = sin(obj.delta);

            sin2beta_num = (sin_e_p^2 + (2 * sin_e_p * sin_e_q * cos(obj.angleMat)) + sin_e_q^2); 
            sin2beta = sin2beta_num ./ ((sin(obj.angleMat)) .^ 2); 
            sin_beta = sqrt(sin2beta); % Np x Nq

            % Pairwise combinations of Rp and context.q 
            y_num = permute(sin_e_p * reshape(obj.q', 1,3,[]) + sin_e_q * Rp, [1,3,2]);
            y_den = sin_beta .* sin(obj.angleMat);

            y = y_num ./ y_den;

            x = obj.pairwiseCross(Rp, obj.q);
            z = cross(x, y);

            cos_beta = cos(asin(sin_beta));

            normals1 = sin_beta .* z + cos_beta .* x;
            normals2 = sin_beta .* z - cos_beta .* x;
        end
        
        function [block] = updateLowerBound(obj, block, thres_stop) 
            %UPDATELOWERBOUND Update block lower bound at stopping threshold
            assert(~isempty(obj.n1_LB) & ~isempty(obj.n2_LB), 'Context was not set');
            
            % Pass near-epipole check option to T search, only at R stopping threshold
            st = StereoT(obj.p, obj.q, obj.n1_LB, obj.n2_LB, obj.t_long_lat, obj.t_half_len, obj.t_half_len_stop, obj.epipole_threshold);
            fprintf("{\n");
            st = st.findSolutions(true); % early_stop = true
            fprintf("}\n");
            block.LB = st.e_max;
            block.patches = st.solutions;
        end
        
        function [block] = updateUpperBound(obj, block, thres_stop) 
            %UPDATEUPPERBOUND Update block upper bound at threshold
            if block.thres > thres_stop
                assert(~isempty(obj.n1_UB) & ~isempty(obj.n2_UB), 'Context was not set');
                st = StereoT(obj.p, obj.q, obj.n1_UB, obj.n2_UB, obj.t_long_lat, obj.t_half_len, obj.t_half_len_stop, -1);
                fprintf("{\n");
                st = st.findSolutions(true); % early_stop = true
                fprintf("}\n");
                block.UB = st.e_max;
            else
                block.UB = block.LB;
            end
        end
        
        function [obj, solutions] = findSolutions(obj)
            init_blk = cube(obj.R_centre, obj.R_sigma);
            obj = obj.bnb(init_blk.subdivide(), obj.thres_stop_R, false); 
            solutions = obj.solutions;
        end
        
    end
    
    methods (Static)
        function [out] = pairwiseCross(v1, v2)
        %PAIRWISECROSS Returns pairwise cross products of two arrays of
        %vectors. Arrays are (num_vectors x dimensions)
            n1 = size(v1, 1);
            n2 = size(v2, 1);
            dim = size(v1,2);
            out = zeros(n1, n2, dim);

            for row = 1 : n1
                out(row,:,:) = cross(repmat(v1(row, :), n2, 1), v2);
            end
        end
        
    end
end

