classdef tTriangle
    %TTriangle Spherical triangle
    properties
        p1       % triangle vertices in counter-clockwise order
        p2
        p3
        
        n12      % normals defining triangle edges (pointing 'into' triangle)
        n23
        n31
        
        centre   
        thres
        sigma
        
        UB       % block upper bound
        LB       % block lower bound
        
        edges_LB    % edges matrix of possible matches at stopping threshold/LB
        edges_UB    % edges matrix of possible matches at threshold/UB
        parent_edges_UB % Edge matrix marking inlier wedges of parent block
    end
    
    methods
        
        function obj = tTriangle(p1, p2, p3, p_edges_UB)
            %TPATCH Construct an instance of this class
            if nargin ~= 0
                obj.p1 = p1 ./ sqrt(sum(p1.^2));
                obj.p2 = p2 ./ sqrt(sum(p2.^2));
                obj.p3 = p3 ./ sqrt(sum(p3.^2));
            end
            obj.n12 = cross(obj.p1, obj.p2);
            obj.n12 = obj.n12 ./ sqrt(sum(obj.n12.^2));
            obj.n23 = cross(obj.p2, obj.p3);
            obj.n23 = obj.n23 ./ sqrt(sum(obj.n23.^2));
            obj.n31 = cross(obj.p3, obj.p1);
            obj.n31 = obj.n31 ./ sqrt(sum(obj.n31.^2));
            
            obj.centre = (p1 + p2 + p3);
            obj.centre = obj.centre ./ sqrt(sum(obj.centre.^2));
            
            obj.sigma = StereoInterface.angles(p1, p2) / 2;
            obj.thres = obj.sigma;
            if exist('p_edges_UB','var') && ~isempty(p_edges_UB)
                obj.parent_edges_UB = p_edges_UB;
            end
        end
        
        function [subblocks] = subdivide(obj)
            %SUBDIVIDE Subdivide triangle and return 4 triangles
            
            % Find midpoints of each triangle edge arc
            m12 = obj.arcMidpoint(obj.p1, obj.p2);
            m23 = obj.arcMidpoint(obj.p2, obj.p3);
            m31 = obj.arcMidpoint(obj.p3, obj.p1);
            
            subblocks = [tTriangle(obj.p1, m12, m31, obj.edges_UB), ...
                         tTriangle(m12, obj.p2, m23, obj.edges_UB), ...
                         tTriangle(m12, m23, m31,    obj.edges_UB), ...
                         tTriangle(m31, m23, obj.p3, obj.edges_UB)];
        end
        
    end
    methods (Static)
        function m = arcMidpoint(p1, p2)
            m = (p1 + p2);
            m = m ./ sqrt(sum(m.^2));
        end
    end
end
