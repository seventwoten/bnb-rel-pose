% Test for t-only case 
% Search with knowledge of ground truth R
% 30 3D points seen in both views, with quantisation noise and no outliers

close all; clear all;

% Set scene variables
load('scenes\s30.mat');

for scene_num = 1:10
scene = s30(scene_num);
diary(['diary\diary_t_tri_',num2str(scene.N),'.',num2str(scene_num),'_',datestr(now,'yyyy-mm-dd','local'),'_',datestr(now,'hh.MM.ss','local'),'.txt'])
diary on;

p = scene.view2;
q = scene.view1;


% Ground truth R and T in required input formats
axis_angle = scene.cam2_aa;
t_long_lat = scene.cam2_tp; 


% Set threshold variables
delta      = 0;             % minimum angular error in Rp and q 
thres_stop = 1/512 * pi;     % stop at this patch size
R_half_len = 1/256 * pi;     % find T given an R block of this size
t_half_len = 1/4 * pi;      % search within an initial T patch of this size 

% Do T search for a particular R block (centred at ground truth)
R_block = RCube(axis_angle, R_half_len);

% Get spherical wedges corresponding to Rp-q pairs
stRT = StereoRT(p, q, [], [], [], [], delta, [], []);
[n1, n2] = stRT.getWedges(R_block, p, q, R_block.thres);

t_list = [tTriangle([0,0,-1],[1,0,0], [0,-1,0]), ...
          tTriangle([-1,0,0],[0,0,-1],[0,-1,0]), ...
          tTriangle([0,0,1],[-1,0,0], [0,-1,0]), ...
          tTriangle([1,0,0],[0,0,1],[0,-1,0]), ...
          tTriangle([1,0,0],[0,0,-1],[0,1,0]), ...
          tTriangle([0,0,-1],[-1,0,0],[0,1,0]), ...
          tTriangle([-1,0,0],[0,0,1],[0,1,0]), ...
          tTriangle([0,0,1],[1,0,0],[0,1,0])];

tic;
stT = StereoT(p, q, n1, n2, t_list, thres_stop, []);
[stT, solutions] = stT.findSolutions(true, false); % For triangles, parallel_mode = false
toc;

fprintf("Number of solutions: %d\n", size(solutions,2));
fprintf("Displaying up to 5:\n");
for s = 1: min(size(solutions,2), 5)
    t_err = StereoInterface.angles(solutions(s).centre, scene.cam2_xyz);
    fprintf("Solution %d:  xyz = [%f %f %f], score: %d, t-error: %f (%f deg)\n", ...
            s, solutions(s).centre, solutions(s).LB, t_err, rad2deg(t_err));
    
    % Plot solution matrices
    figure, imagesc(solutions(s).edges_LB);
end

fprintf("Ground truth: aa = [%f %f %f], rpy = [%f %f %f], theta-phi = [%f %f]\n", ...
        scene.cam2_aa, scene.cam2_rpy, scene.cam2_xyz);

diary off;

end
