classdef tPatchList < tPatch
    %TPATCHLIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nb         % number of blocks represented in list
    end
    
    methods
        function obj = tPatchList(c, s, R, p_edges_UB)
            %TPATCHLIST Construct an instance of this class
            %   c is a matrix formed by concatenating patches centres (Nb x 2)
            %   s is the patch half-length shared by all patches
            if nargin == 0
                c = [];
                s = [];
            end
            if ~exist('R', 'var')
                R = [];
            end
            if ~exist('p_edges_UB', 'var')
                p_edges_UB = [];
            end
            obj = obj@tPatch(c, s, R, p_edges_UB);
            obj.Nb = size(c, 1);
            obj.LB = zeros(obj.Nb, 1);
            obj.UB = zeros(obj.Nb, 1);
        end
        
        function [subblocks] = subdivide(obj, mask)
            %SUBDIVIDE Return new tPatchList containing subblocks of blocks 
            %indicated by 1 in input mask (Nb x 1)
            if ~exist('mask','var')
                mask = true(obj.Nb, 1);
            end
            
            if ~isempty(obj.R_tilt)
                obj = obj.rotate(obj.R_tilt');
            end
            shifts = [-1, -1 ;
                      -1,  1 ;
                       1, -1 ;
                       1,  1 ]' ;
            
            sigma_new  = obj.sigma/2;
            
            % Select centres to subdivide where mask == true
            % Transpose centres to column vectors for easier broadcasting
            parent_blocks = obj.centre(mask, :)';
            parent_blocks = reshape(parent_blocks, 2,1,[]);
            
            centres_new = parent_blocks + sigma_new * shifts;
            centres_new = reshape(centres_new, 2, [])';

            parent_edges_UB = obj.edges_UB(:, :, mask);
            parent_edges_UB = repelem(parent_edges_UB, 1, 1, 4);
            
            subblocks = tPatchList(centres_new, sigma_new, obj.R_tilt, parent_edges_UB);
            
        end
        
        function obj = removePatches(obj, mask) 
            remaining = ~mask;
            obj.centre = obj.centre(remaining, :);
            obj.centre_xyz = obj.centre_xyz(remaining, :);
            obj.UB = obj.UB(remaining, :);
            obj.LB = obj.LB(remaining, :);
            obj.thres = obj.thres(remaining, :);
            obj.Nb = size(obj.thres, 1);
        end
    end
end

