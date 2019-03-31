function [result] = angles(Rp, q)
%MULTIPLY Finds pairwise angles between two arrays of 3D vectors
%   Detailed explanation goes here

Rp_mag = sqrt(sum(Rp.^2,2)); % N x 1
q_mag  = sqrt(sum( q.^2,2)); % N x 1
divisor = Rp_mag * q_mag';   % N x N

cos_angle = (Rp * q') ./ divisor;
result = abs(acos(cos_angle));

end


