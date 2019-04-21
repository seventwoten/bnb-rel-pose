classdef cubeRT < cube
    %cubeRT 3D block for searching R and t
    properties
        n1
        n2
    end
    
    methods
        function obj = cubeRT(c, s, p)
            %CUBERT Construct an instance of this class
            if nargin == 0
                c = [];
                s = [];
                p = [];
            end
            obj = obj@cube(c, s, p);
        end
        
        function [obj] = setContext(obj, context, delta)
            %SETCONTEXT Compute (n1, n2) from (R, p, q, angles, thres, delta)
            R  = aa2mat(obj.centre);
            Rp = (R * context.p')';
            obj.angleMat = angles(Rp, context.q);
            
            sin_e_p = sin(delta + obj.thres);
            sin_e_q = sin(delta);

            sin2beta_num = (sin_e_p^2 + (2 * sin_e_p * sin_e_q * cos(obj.angleMat)) + sin_e_q^2); 
            sin2beta = sin2beta_num ./ ((sin(obj.angleMat)) .^ 2); 
            sin_beta = sqrt(sin2beta); % Np x Nq

            % Pairwise combinations of Rp and context.q 
            y_num = permute(sin_e_p * reshape(context.q', 1,3,[]) + sin_e_q * Rp, [1,3,2]);
            y_den = sin_beta .* sin(obj.angleMat);

            y = y_num ./ y_den;

            x = obj.pairwiseCross(Rp, context.q);
            z = cross(y, x);

            cos_beta = cos(asin(sin_beta));

            obj.n1 = sin_beta .* z + cos_beta .* x;
            obj.n2 = sin_beta .* z - cos_beta .* x;
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

