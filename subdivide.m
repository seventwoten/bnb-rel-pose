function [centre_new, sigma_new] = subdivide(centre, sigma)
%SUBDIVIDE Summary of this function goes here
%   Detailed explanation goes here
shifts = [-1, -1, -1 ;
          -1, -1,  1 ;
          -1,  1, -1 ;
          -1,  1,  1 ;
           1, -1, -1 ;
           1, -1,  1 ;
           1,  1, -1 ;
           1,  1,  1  ];
       
sigma_new  = sigma/2;
centre_new = centre + sigma_new * shifts;

end
