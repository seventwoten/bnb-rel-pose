classdef StereoRT < StereoInterface
    %STEREORT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p_known
        q_known
        possibleMatches    % Np x Nq mask for including only possible matches
        
        R_list             % List of RCubes
        thres_stop_R
        
        t_list             % Default list of tPatches
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
        function obj = StereoRT(p, q, R_list, R_sigma_stop, t_list, t_half_len_stop, delta, epipole_threshold, corr_indices, possible_matches)
            %STEREORT Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@StereoInterface(p, q);
            obj.thres_stop_R = sqrt(3) * R_sigma_stop;
            obj.t_half_len_stop = t_half_len_stop;
            obj.delta = delta;
            obj.epipole_threshold = epipole_threshold;
            
            if exist('possible_matches', 'var') && ~isempty(possible_matches)
                assert(size(possible_matches, 1) == obj.Np && size(possible_matches, 2) == obj.Nq, 'PossibleMatches is the wrong size');
                obj.possibleMatches = possible_matches;
            end
            
            % Separate out known correspondence
            if ~isempty(corr_indices)
                obj.p_known = p(corr_indices(:,1), :);
                obj.q_known = q(corr_indices(:,2), :);
                obj.p(corr_indices(:,1), :) = [];
                obj.q(corr_indices(:,2), :) = [];
                obj.Np = size(obj.p, 1);
                obj.Nq = size(obj.q, 1);
                if ~isempty(obj.possibleMatches)
                    obj.possibleMatches(corr_indices(:,1), :) = [];
                    obj.possibleMatches(:, corr_indices(:,2)) = [];
                end
            end
            
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
            [obj.n1_LB, obj.n2_LB] = getWedges(obj, block, obj.p, obj.q, obj.thres_stop_R);
        end
        
        function [normals1, normals2] = getWedges(obj, block, p, q, threshold_R)
            R  = block.aa2mat();
            Rp = (R * p')';
            obj.angleMat = obj.angles(Rp, q);
            
            e_p = obj.delta + threshold_R;
            e_q = obj.delta;
            
            sin_e_p = sin(e_p);
            sin_e_q = sin(e_q);

            sq_terms = sin_e_p .* sin_e_p + sin_e_q .* sin_e_q;
            coeff_cos = 2 .* sin_e_p .* sin_e_q;
            sin2beta_num = sq_terms + (coeff_cos .* cos(obj.angleMat));
            sin2beta = sin2beta_num ./ ((sin(obj.angleMat)) .^ 2); 
            sin_beta = sqrt(sin2beta); % Np x Nq

            % Pairwise combinations of Rp and q
            y_num = permute(sin_e_p * reshape(q', 1,3,[]) + sin_e_q * Rp, [1,3,2]);
            y_den = sin_beta .* sin(obj.angleMat);

            y = y_num ./ y_den;

            x = obj.pairwiseCross(Rp, q);
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
        
        function [t_list] = wedge2patches(obj, n1_wedge, n2_wedge, min_size)
            %SETCONTEXT Convert a single wedge into a list of tPatches.
            %   Align wedge to z-axis, divide into patches, and rotate back.
            
            % Normalise and reshape wedge normals to unit row vectors
            n1_wedge = reshape(n1_wedge, [], 3) ./ sqrt(sum(n1_wedge.^2));
            n2_wedge = reshape(n2_wedge, [], 3) ./ sqrt(sum(n2_wedge.^2));
            
            % Find wedge centre and corner (unit vectors), and wedge width
            wedge_centre = (n1_wedge + n2_wedge) / 2; 
            wedge_centre = wedge_centre ./ sqrt(sum(wedge_centre.^2));
            
            wedge_corner = cross(wedge_centre, n1_wedge); 
            wedge_corner = wedge_corner ./ sqrt(sum(wedge_corner.^2));
            
            wedge_width  = abs(acos(n1_wedge * -n2_wedge'));
            
            % Find rotation to align wedge corner to z-axis
            v   = cross(wedge_corner, [0,0,1]);
            v_x = [   0, -v(3),  v(2); 
                   v(3),     0, -v(1); 
                  -v(2),  v(1),     0];

            c = wedge_corner * [0,0,1]';
            if abs(c-1) < 1e-8 || abs(c+1) < 1e-8 % No rotation needed
                R_align = eye(3);
            else
                R_align = eye(3) + v_x + (v_x * v_x ./ (1+c));
            end
            
            % Find wedge_centre coords after aligning wedge corner to z-axis
            wedge_centre_aln = R_align * wedge_centre';
            
            % Find number and size of patches along wedge width
            if ~exist('min_size', 'var') || isempty(min_size)
                min_size = max(pi/16, obj.t_half_len_stop);
            end
            num_theta      = ceil( wedge_width / (2*min_size) );
            patch_half_len = min_size;
            while mod(num_theta, 2) == 0
                num_theta      = num_theta  / 2;
                patch_half_len = patch_half_len * 2;
            end
            
            % Find number of patches along wedge length
            num_phi = ceil( pi / (2*patch_half_len) );
            
            % Initialise tPatches aligned to z-axis, then rotate back
            rtp = tPatch.cartesian2Spherical(wedge_centre_aln(1),wedge_centre_aln(2),wedge_centre_aln(3)); 
            thetas = rtp(2) + (-floor(num_theta/2):floor(num_theta/2)) * patch_half_len * 2;
            phis   = rtp(3) + (-floor(num_phi/2)  :floor(num_phi/2))   * patch_half_len * 2;
            
            [t,p] = meshgrid(thetas,phis);
            t_centres = [t(:), p(:)];
            t_list = tPatchList(t_centres, patch_half_len, R_align');
        end
        
        function t_list = updateTList(obj, block, t_list, thres_R)
            % Default case: Return given t_list
            
            % Single known correspondence: Find known p-q wedge, use as new t_list
            if ~isempty(obj.p_known) && ~isempty(obj.q_known)
                [n1_wedge, n2_wedge] = obj.getWedges(block, obj.p_known, obj.q_known, thres_R);
                
                % If wedge normals are valid, replace t_list with wedge patches
                if ~all(n1_wedge == 0) && ~all(n2_wedge == 0)
                    t_list = obj.wedge2patches(n1_wedge, n2_wedge);
                end
            end
        end
        
        function [block] = updateLowerBound(obj, block, thres_stop) 
            %UPDATELOWERBOUND Update block lower bound at stopping threshold
            assert(~isempty(obj.n1_LB) & ~isempty(obj.n2_LB), 'Context was not set');
            
            % Update T search list if known correspondence exists
            tlist = obj.updateTList(block, obj.t_list, obj.thres_stop_R);
            
            % Pass near-epipole check option to T search, only at R stopping threshold
            st = StereoT(obj.p, obj.q, obj.n1_LB, obj.n2_LB, tlist, obj.t_half_len_stop, obj.epipole_threshold, obj.possibleMatches);
            fprintf("{\n");
            st = st.findSolutions(true); % early_stop = true
            fprintf("}\n");
            block.LB = st.e_max;
            block.patches = st.solutions;
        end
        
        function [block] = updateUpperBound(obj, block)
            %UPDATEUPPERBOUND Update block upper bound at threshold
            [obj.n1_UB, obj.n2_UB] = getWedges(obj, block, obj.p, obj.q, block.thres);
            
            % Update T search list if known correspondence exists
            tlist = obj.updateTList(block, obj.t_list, block.thres);
            
            st = StereoT(obj.p, obj.q, obj.n1_UB, obj.n2_UB, tlist, obj.t_half_len_stop, -1, obj.possibleMatches);
            fprintf("{\n");
            st = st.findSolutions(true); % early_stop = true
            fprintf("}\n");
            block.UB = st.e_max;
        end
        
        function [obj, solutions] = findSolutions(obj, early_stop)
            % Default: Return all surviving blocks as solutions
            if ~exist('early_stop', 'var') || isempty(early_stop)
                early_stop = false;
            end
            
            init_list = obj.R_list;
            if numel(init_list) == 1
                init_list = init_list.subdivide();
            end
            obj = obj.bnb(init_list, obj.thres_stop_R, early_stop);
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

