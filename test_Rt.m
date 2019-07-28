% Test for unknown R and t, with up to 1 known correspondence
% 30 3D points seen in both views, with quantisation noise and no outliers

% Ground truth: 
% rpy: [0.2,-0.1,0.1] -> axis-angle as 3d point in R ball: [-0.1095839, 0.089599, 0.2046666]
% xyz: [-0.3,-0.4, 0] -> theta-phi as longitude-latitude: [-2.2142974, 1.5707963]

close all; clear all;

load('scenes\s30.mat');
scene = s30(1);

p = scene.view2;
q = scene.view1;
known_corr = [1,1];    % set to [] for no known correspondences
diary_corr = '1c_';    % set to [] to turn off diary name

% Set input variables
R_list       = [RCube([-3/8 * pi, -3/8 * pi, 0], pi/8), RCube([    -pi/8, -3/8 * pi, 0], pi/8), ...
                RCube([-3/8 * pi,     -pi/8, 0], pi/8), RCube([    -pi/8,     -pi/8, 0], pi/8), ...
                RCube([-3/8 * pi,      pi/8, 0], pi/8), RCube([    -pi/8,      pi/8, 0], pi/8), ...
                RCube([-3/8 * pi,  3/8 * pi, 0], pi/8), RCube([    -pi/8,  3/8 * pi, 0], pi/8), ...
                RCube([     pi/8, -3/8 * pi, 0], pi/8), RCube([ 3/8 * pi, -3/8 * pi, 0], pi/8), ...
                RCube([     pi/8,     -pi/8, 0], pi/8), RCube([ 3/8 * pi,     -pi/8, 0], pi/8), ...
                RCube([     pi/8,      pi/8, 0], pi/8), RCube([ 3/8 * pi,      pi/8, 0], pi/8), ...
                RCube([     pi/8,  3/8 * pi, 0], pi/8), RCube([ 3/8 * pi,  3/8 * pi, 0], pi/8)];

t_list       = [tPatch([0,pi/2], pi/2), tPatch([pi,pi/2], pi/2)];

delta        = 0;             % minimum angular error in Rp and q 
thres_stop_R = 1/32 * pi;     % Stop when cube diagonal drops below this value
thres_stop_t = 1/64 * pi;     % Stop when patch diagonal drops below this value
epipole_thres = 0.7;          % Reject points that match more than this fraction of points


diary([['diary\diary_Rt_', diary_corr],datestr(now,'yyyy-mm-dd','local'),'_',datestr(now,'hh.MM.ss','local'),'.txt'])
diary on;

fprintf("test_Rt: \n");
fprintf("delta         =  %f \n",   delta);
fprintf("R_half_len    = %s pi \n", rats(R_list(1).sigma/pi,5));
fprintf("thres_stop_R  = %s pi \n", rats(thres_stop_R/pi,5));
fprintf("t_half_len    = %s pi \n", rats(t_list(1).sigma/pi,5));
fprintf("thres_stop_t  = %s pi \n", rats(thres_stop_t/pi,5));
fprintf("epipole_thres =  %f \n\n", epipole_thres);

tic;

st = StereoRT(p, q, R_list, thres_stop_R, t_list, thres_stop_t, delta, epipole_thres, known_corr);
[st, solutions] = st.findSolutions();

toc;

fprintf("Number of solutions: %d\n", size(solutions,2));
fprintf("Displaying up to 5:\n");

num_sol = size(solutions,2);
for s = 1: min(num_sol, 5)
    fprintf("Solution %d: [%f %f %f], sigma: %s pi, score: %d rpy: [%f %f %f]\n", ...
            s, solutions(s).centre, rats(solutions(s).sigma/pi,5), solutions(s).LB, ...
            RCube.R2rpy(solutions(s).aa2mat()));
    % For each R, plot 1 solution matrix representing one possible t
    figure, imagesc(solutions(s).patches(1).edges_stop);
end

fprintf("Ground truth: axis-angle = [%f %f %f], rpy = [%f %f %f], theta-phi = [%f %f]\n", ... 
        scene.cam2_aa, scene.cam2_rpy, scene.cam2_tp);

diary off;
