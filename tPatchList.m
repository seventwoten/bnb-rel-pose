classdef tPatchList < tPatch
    %TPATCHLIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nb         % number of blocks represented in list
    end
    
    methods
        function obj = tPatchList(c, s)
            %TPATCHLIST Construct an instance of this class
            %   c is a matrix formed by concatenating patches centres (Nb x 2)
            %   s is the patch half-length shared by all patches
            if nargin == 0
                c = [];
                s = [];
            end
            obj = obj@tPatch(c, s);
            obj.Nb = size(c, 1);
            obj.LB = zeros(obj.Nb, 1);
            obj.UB = zeros(obj.Nb, 1);
            
        end
        
        function [subblocks] = subdivide(obj, mask)
            %SUBDIVIDE Return new tPatchList containing subblocks of blocks 
            %indicated by 1 in input mask (Nb x 1)
            if nargin == 1
                mask = true(obj.Nb, 1);
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

            subblocks = tPatchList(centres_new, sigma_new);
            
        end
    end
end

