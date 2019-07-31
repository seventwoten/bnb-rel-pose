% Test for unknown R and t, with up to 1 known correspondence
% 30 3D points seen in both views, with quantisation noise and no outliers

close all; clear all;

% Set scene variables
load('scenes\s30.mat');
scene_num = 1;
scene = s30(scene_num);

p = scene.view2;
q = scene.view1;

known_corr = [1,1];    % set to [] for no known correspondences
if (size(known_corr, 1) == 1); diary_corr = '1c_'; else; diary_corr = []; end

% Set search variables
R_list       = [RCube([-3/8 * pi, -3/8 * pi, 0], pi/8), RCube([    -pi/8, -3/8 * pi, 0], pi/8), ...
                RCube([-3/8 * pi,     -pi/8, 0], pi/8), RCube([    -pi/8,     -pi/8, 0], pi/8), ...
                RCube([-3/8 * pi,      pi/8, 0], pi/8), RCube([    -pi/8,      pi/8, 0], pi/8), ...
                RCube([-3/8 * pi,  3/8 * pi, 0], pi/8), RCube([    -pi/8,  3/8 * pi, 0], pi/8), ...
                RCube([     pi/8, -3/8 * pi, 0], pi/8), RCube([ 3/8 * pi, -3/8 * pi, 0], pi/8), ...
                RCube([     pi/8,     -pi/8, 0], pi/8), RCube([ 3/8 * pi,     -pi/8, 0], pi/8), ...
                RCube([     pi/8,      pi/8, 0], pi/8), RCube([ 3/8 * pi,      pi/8, 0], pi/8), ...
                RCube([     pi/8,  3/8 * pi, 0], pi/8), RCube([ 3/8 * pi,  3/8 * pi, 0], pi/8)];
t_list       = [tPatch([0,pi/2], pi/2), tPatch([pi,pi/2], pi/2)];

delta         = 0;             % minimum angular error in Rp and q 
thres_stop_R  = 1/256 * pi;   % Stop when cube diagonal drops below this value
thres_stop_t  = 1/256 * pi;   % Stop when patch diagonal drops below this value
epipole_thres = 0.7;          % Reject points that match more than this fraction of points
early_stop    = true;         % Set to true to return first viable solution

% Print output to diary file
timestamp = [datestr(now,'yyyy-mm-dd','local'),'_',datestr(now,'hh.MM.ss','local')];
diary([['diary\diary_Rt_', num2str(scene.N),'.',num2str(scene_num),'_', diary_corr], timestamp, '.txt'])
diary on;

fprintf("test_Rt: \n");
fprintf("delta         =  %f \n",   delta);
fprintf("R_half_len    = %s pi \n", rats(R_list(1).sigma/pi,5));
fprintf("thres_stop_R  = %s pi \n", rats(thres_stop_R/pi,5));
fprintf("t_half_len    = %s pi \n", rats(t_list(1).sigma/pi,5));
fprintf("thres_stop_t  = %s pi \n", rats(thres_stop_t/pi,5));
fprintf("epipole_thres =  %f \n",   epipole_thres);
fprintf("early_stop    =  %d \n\n", early_stop);


tic;

st = StereoRT(p, q, R_list, thres_stop_R, t_list, thres_stop_t, delta, epipole_thres, known_corr);
[st, solutions] = st.findSolutions(early_stop);
fprintf("--------");

toc;

fprintf("Number of solutions: %d (Displaying up to 5)\n", size(solutions,2));

num_sol = size(solutions,2);
for s = 1: min(num_sol, 5)
    % Report R and t angular error
    R_err = StereoInterface.angles([0,0,1]*solutions(s).aa2mat()', [0,0,1] * scene.cam2_R');
    t_err = StereoInterface.angles(solutions(s).patches(1).centre_xyz, scene.cam2_xyz);
    fprintf("Solution %d:   aa = [%f %f %f], rpy: [%f %f %f], theta-phi = [%f %f], score: %d, R-error: %f (%f deg), t-error: %f (%f deg)\n", ...
            s, solutions(s).centre, RCube.R2rpy(solutions(s).aa2mat()), solutions(s).patches(1).centre, solutions(s).LB, ...
            R_err, rad2deg(R_err), t_err, rad2deg(t_err));
    % For each R, plot 1 solution matrix representing one possible t
    figure, imagesc(solutions(s).patches(1).edges_stop);
end

fprintf("Ground truth: aa = [%f %f %f], rpy = [%f %f %f], theta-phi = [%f %f]\n", ...
        scene.cam2_aa, scene.cam2_rpy, scene.cam2_tp);

diary off;

% Save solutions to .mat file
save(['solutions\Rt_', num2str(scene.N),'.',num2str(scene_num),'_',diary_corr, timestamp, '.mat'],'solutions');
