% Test for R-only case
% 30 3D points seen in both views, with quantisation noise and no outliers

% Ground truth: 
% R  = [ 0.97517033 -0.20744485  0.07752075;
%        0.19767681  0.97319025  0.11757815;
%       -0.09983342 -0.09933467  0.99003329]

% roll_pitch_yaw = [0.2, -0.1, 0.1]

close all; clear all;
diary(['diary\diary_R_',datestr(now,'yyyy-mm-dd','local'),'_',datestr(now,'hh.MM.ss','local'),'.txt']);
diary on;

p = [ 0.0548, -0.2238,  1   ;
      0.0817,  0.2685,  1   ;
      0.4416,  0.2896,  1   ;
      0.207 ,  0.1147,  1   ;
      0.26  , -0.1423,  1   ;
      0.2682, -0.3267,  1   ;
      0.5068, -0.5963,  1   ;
     -0.298 ,  0.4571,  1   ;
      0.0402, -0.1252,  1   ;
     -0.3346, -0.0481,  1   ;
      0.2475, -0.1129,  1   ;
      0.2882,  0.3811,  1   ;
      0.3109,  0.014 ,  1   ;
      0.4682,  0.5085,  1   ;
      0.3719,  0.5554,  1   ;
     -0.2343, -0.2506,  1   ;
      0.4168,  0.1615,  1   ;
      0.5082,  0.1688,  1   ;
     -0.0134, -0.2721,  1   ;
     -0.4133, -0.1682,  1   ;
      0.6528, -0.1418,  1   ;
      0.372 ,  0.225 ,  1   ;
      0.1272,  0.2764,  1   ;
     -0.2483, -0.2032,  1   ;
     -0.1404, -0.2542,  1   ;
      0.3851,  0.2668,  1   ;
      0.4205,  0.2952,  1   ;
     -0.3189,  0.0342,  1   ;
      0.1244,  0.0579,  1   ;
      0.3902, -0.3613,  1   ];


q = [ 0.0196, -0.2977,  1   ;
     -0.0386,  0.2204,  1   ;
      0.3546,  0.3522,  1   ;
      0.1467,  0.122 ,  1   ;
      0.2519, -0.1245,  1   ;
      0.2178, -0.4023,  1   ;
      0.4177, -0.7346,  1   ;
     -0.469 ,  0.3127,  1   ;
      0.0346, -0.1407,  1   ;
     -0.3649, -0.1725,  1   ;
      0.2463, -0.081 ,  1   ;
      0.1479,  0.3842,  1   ;
      0.2904,  0.0623,  1   ;
      0.279 ,  0.5216,  1   ;
      0.1462,  0.5272,  1   ;
     -0.262 , -0.3886,  1   ;
      0.3309,  0.1799,  1   ;
      0.4196,  0.1962,  1   ;
     -0.0035, -0.3102,  1   ;
     -0.5476, -0.4823,  1   ;
      0.5175, -0.2502,  1   ;
      0.3155,  0.2973,  1   ;
      0.0398,  0.2842,  1   ;
     -0.229 , -0.2725,  1   ;
     -0.1275, -0.3143,  1   ;
      0.2023,  0.1866,  1   ;
      0.2851,  0.2926,  1   ;
     -0.3378, -0.0493,  1   ;
      0.0581,  0.0255,  1   ;
      0.4144, -0.3264,  1   ];

tic;

stR = StereoR(p, q, [0,0,0], pi, 1/32 * pi, 1/128 * pi);
[stR, solutions] = stR.findSolutions();

toc;

for s = 1: size(solutions,2)
    fprintf("Solution %d: [%d %d %d], sigma: %d pi, score: %d\n R:", s, solutions(s).centre, solutions(s).sigma/pi, solutions(s).LB);
    R = solutions(s).aa2mat();
    disp(R);
    fprintf("rpy: [%d %d %d]\n\n", R2rpy(R));
    
    % Plot solution matrices
    figure, imagesc(solutions(s).edges_stop);
end

diary off;
