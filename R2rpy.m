function [rpy] = R2rpy(R)
%R2RPY Summary of this function goes here
%   Detailed explanation goes here
    
    sy = sqrt(R(1,1) * R(1,1) + R(2,1) * R(2,1));
    singular = sy < 1e-6;
    
    if ~singular
        r = atan2( R(2,1), R(1,1));
        p = atan2( R(3,2), R(3,3));
        y = atan2(-R(3,1), sy);
    else
        p = atan2(-R(2,3), R(2,2));
        y = atan2(-R(3,1), sy);
        r = 0;
    end
        
    rpy = [r, p, y];
end
