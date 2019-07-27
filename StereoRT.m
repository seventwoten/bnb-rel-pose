classdef StereoRT < StereoInterface
    %STEREORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        R_list             % List of RCubes
        thres_stop_R
        
        t_list             % List of tPatches
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
        function obj = StereoRT(p, q, R_list, R_sigma_stop, t_list, t_half_len_stop, delta, epipole_threshold)
            %STEREORT Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@StereoInterface(p, q);
            obj.thres_stop_R = sqrt(3) * R_sigma_stop;
            obj.t_half_len_stop = t_half_len_stop;
            obj.delta = delta;
            obj.epipole_threshold = epipole_threshold;
            
            if isempty(R_list)
                % default R range to search
                obj.R_list = RCube([0, 0, 0], pi);
            else
                obj.R_list = R_list;
            end
            
            if isempty(t_list)
                % default t range to search
                obj.t_list = [tPatch([0, 0], pi/2), tPatch([0, pi], pi/2)];
            else
                obj.t_list = t_list;
            end
            
        end
        
        function [obj] = setContext(obj, block)
            %SETCONTEXT Compute lower bound (n1, n2) from (R, p, q, angles, thres, delta)
            [obj.n1_LB, obj.n2_LB] = getWedges(obj, block, obj.thres_stop_R);
        end
        
        function [normals1, normals2] = getWedges(obj, block, threshold_R)
            R  = block.aa2mat();
            Rp = (R * obj.p')';
            obj.angleMat = obj.angles(Rp, obj.q);
            
            e_p = obj.delta + threshold_R;
            e_q = obj.delta;
            
            sin_e_p = sin(e_p);
            sin_e_q = sin(e_q);

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
            
            % Zero out normals where Rp and q circular ranges overlap 
            % (wedge constraint does not hold)
            has_overlap = obj.angleMat < e_p + e_q;
            normals1(has_overlap(:,:,[1,1,1])) = 0;
            normals2(has_overlap(:,:,[1,1,1])) = 0;
            
        end
        
        function [block] = updateLowerBound(obj, block, thres_stop) 
            %UPDATELOWERBOUND Update block lower bound at stopping threshold
            assert(~isempty(obj.n1_LB) & ~isempty(obj.n2_LB), 'Context was not set');
            
            % Pass near-epipole check option to T search, only at R stopping threshold
            st = StereoT(obj.p, obj.q, obj.n1_LB, obj.n2_LB, obj.t_list, obj.t_half_len_stop, obj.epipole_threshold);
            fprintf("{\n");
            st = st.findSolutions(true); % early_stop = true
            fprintf("}\n");
            block.LB = st.e_max;
            block.patches = st.solutions;
        end
        
        function [block] = updateUpperBound(obj, block)
            %UPDATEUPPERBOUND Update block upper bound at threshold
            [obj.n1_UB, obj.n2_UB] = getWedges(obj, block, block.thres);
            st = StereoT(obj.p, obj.q, obj.n1_UB, obj.n2_UB, obj.t_list, obj.t_half_len_stop, -1);
            fprintf("{\n");
            st = st.findSolutions(true); % early_stop = true
            fprintf("}\n");
            block.UB = st.e_max;
        end
        
        function [obj, solutions] = findSolutions(obj)
            init_list = obj.R_list;
            if numel(init_list) == 1
                init_list = init_list.subdivide();
            end
            obj = obj.bnb(init_list, obj.thres_stop_R, false); 
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

