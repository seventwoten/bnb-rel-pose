% Test for t-only case 
% Search is close to ground truth R and T
% 30 3D points seen in both views, with quantisation noise and no outliers

% Ground truth: 
% rpy: [0.2,-0.1,0.1] -> axis-angle as 3d point in R ball: [-0.1095839, 0.089599, 0.2046666]
% xyz: [-0.3,-0.4, 0] -> theta-phi as longitude-latitude: [-2.2142974, 1.5707963]

close all; clear all;
diary(['diary\diary_t_',datestr(now,'yyyy-mm-dd','local'),'_',datestr(now,'hh.MM.ss','local'),'.txt'])
diary on;

p = [-0.0596,  0.2702, 1.  ;
      0.3342, -0.1944, 1.  ;
      0.242 , -0.1989, 1.  ;
     -0.1747, -0.2492, 1.  ;
      0.3273,  0.3383, 1.  ;
     -0.2507,  0.0308, 1.  ;
      0.1016,  0.5406, 1.  ;
     -0.5639, -0.2959, 1.  ;
      0.3975, -0.0002, 1.  ;
     -0.3473, -0.3172, 1.  ;
     -0.3618,  0.2303, 1.  ;
     -0.119 ,  0.2816, 1.  ;
      0.1135, -0.0818, 1.  ;
      0.1882,  0.0011, 1.  ;
     -0.7126, -0.3211, 1.  ;
     -0.0104,  0.3141, 1.  ;
     -0.4788, -0.0284, 1.  ;
      0.8599, -0.6003, 1.  ;
     -0.0263,  0.8015, 1.  ;
      0.184 , -0.0701, 1.  ;
      0.0721,  0.3367, 1.  ;
     -0.1654,  0.44  , 1.  ;
      0.1912, -0.1834, 1.  ;
     -0.1941,  0.0929, 1.  ;
      0.6122, -0.1875, 1.  ;
      0.3611,  0.1535, 1.  ;
     -0.3344, -0.2078, 1.  ;
     -0.3695,  0.2464, 1.  ;
     -0.0065,  0.3912, 1.  ;
      0.3296, -0.1865, 1.  ];

q = [-0.1439,  0.2391, 1. ;
      0.2782, -0.241 , 1. ;
      0.2164, -0.22  , 1. ;
     -0.1804, -0.3419, 1. ;
      0.1984,  0.347 , 1. ;
     -0.2905, -0.0622, 1. ;
     -0.0908,  0.5028, 1. ;
     -0.6058, -0.5603, 1. ;
      0.3073, -0.0367, 1. ;
     -0.3319, -0.4418, 1. ;
     -0.5738, -0.0664, 1. ;
     -0.2083,  0.2329, 1. ;
      0.0681, -0.1252, 1. ;
      0.143 , -0.007 , 1. ;
     -0.7273, -0.6051, 1. ;
     -0.1084,  0.2911, 1. ;
     -0.5479, -0.2438, 1. ;
      0.8283, -0.643 , 1. ;
     -0.3766,  0.6421, 1. ;
      0.1386, -0.0976, 1. ;
     -0.0261,  0.3398, 1. ;
     -0.3094,  0.3628, 1. ;
      0.1729, -0.1999, 1. ;
     -0.2517,  0.0087, 1. ;
      0.472 , -0.3153, 1. ;
      0.3021,  0.1981, 1. ;
     -0.3239, -0.3139, 1. ;
     -0.4886,  0.0763, 1. ;
     -0.1768,  0.3004, 1. ;
      0.3177, -0.1731, 1. ];

% Ground truth R and T in required input formats
axis_angle = [-0.1095839, 0.089599, 0.2046666];
t_long_lat = [-2.2142974, 1.5707963];  

tic;

% Set threshold variables
delta      = 0;             % minimum angular error in Rp and q 
thres_stop = 1/64 * pi;     % stop at this patch size
R_half_len = 1/64 * pi;     % find T given an R block of this size
t_half_len = 1/2 * pi;      % search within an initial T patch of this size

% Do T search for a particular R block (centred at ground truth)
R_block = RCube(axis_angle, R_half_len);

% Get spherical wedges corresponding to Rp-q pairs
stRT = StereoRT(p, q, [], [], [], [], [], [], delta, []);
[n1, n2] = stRT.getWedges(R_block, R_block.thres);

% Start T search within initial patch (centred at ground truth)
stT = StereoT(p, q, n1, n2, t_long_lat, t_half_len, thres_stop, []);
[stT, solutions] = stT.findSolutions(false);

toc;

fprintf("Number of solutions: %d\n", size(solutions,2));
fprintf("Displaying up to 5:\n");
for s = 1: min(size(solutions,2), 5)
    fprintf("Solution %d: [%d %d], sigma: %d pi, score: %d\n", s, solutions(s).centre, solutions(s).sigma/pi, solutions(s).LB);
    % Plot solution matrices
    figure, imagesc(solutions(s).edges_LB);
end

diary off;
