function [R] = aa2mat(b_xyz)
%AA2MAT Given a 3d point in R-space ball, returns R matrix
%   b_xyz is 1 x 3
alpha = sqrt(sum(b_xyz.^2,2));
u = b_xyz / alpha;

c = cos(alpha);
s = sin(alpha);
t = 1-c;

U_outer = u' * u;
U_as = [      0, -u(1,3),  u(1,2) ; 
         u(1,3),       0, -u(1,1) ; 
        -u(1,2),  u(1,1),      0 ];

R = c * eye(3) + t * U_outer + s * U_as;

end

