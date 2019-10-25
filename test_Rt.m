% Test for angel scene with unknown R and t, 2 known correspondences
% >800 points in two views 

close all; clear all;

% Set scene variables p and q (predefined in a script)
% angelpts;
angelpts_s;

% Experiment parameters
known_corr = [701, 685; 122, 123];    % set to [] for no known correspondences
                                      % angelpts:   [777, 762; 136, 140];
                                      % angelpts_s: [701, 685; 122, 123];
possible_matches = []; %blkdiag(true, true, true(scene.N-2)); %set to [] to allow all pairings

expt_name  = [];
if (size(known_corr, 1) > 0); expt_name = [num2str(size(known_corr, 1)), 'r_']; end
if (size(possible_matches, 1) ~= 0); expt_name = '1f_'; end

% Set search variables
R_list       = [RCube([0, -3/8 * pi, 0], pi/8), RCube([0, -pi/8, 0], pi/8); 
                RCube([0,  pi/8, 0], pi/8), RCube([0,  3/8 * pi, 0], pi/8)];
t_list       = [tPatch([0,pi/2], pi/2), tPatch([pi,pi/2], pi/2)];

delta         = 0.00873;      % minimum angular error in Rp and q
thres_stop_R  = 0;            % Stop when cube half-length drops below this value
thres_stop_t  = 1/512 * pi;   % Stop when patch half-length drops below this value
epipole_thres = 0.7;          % Reject points that match more than this fraction of points
early_stop    = true;         % Set to true to return first viable solution

% Print output to diary file
timestamp = [datestr(now,'yyyy-mm-dd','local'),'_',datestr(now,'hh.MM.ss','local')];
diary([['diary\diary_Rt_angel1_', expt_name], timestamp, '.txt'])
diary on;

fprintf("test_Rt: \n");
fprintf("delta         =  %f \n",   delta);
fprintf("R_half_len    = %s pi \n", rats(R_list(1).sigma/pi,7));
fprintf("thres_stop_R  = %s pi \n", rats(thres_stop_R/pi,   7));
fprintf("t_half_len    = %s pi \n", rats(t_list(1).sigma/pi,7));
fprintf("thres_stop_t  = %s pi \n", rats(thres_stop_t/pi,   7));
fprintf("epipole_thres =  %f \n",   epipole_thres);
fprintf("early_stop    =  %d \n\n", early_stop);


tic;

st = StereoRT(p, q, R_list, thres_stop_R, t_list, thres_stop_t, delta, epipole_thres, known_corr, possible_matches);
[st, solutions] = st.findSolutions(early_stop);
fprintf("--------");

toc;

fprintf("Number of solutions: %d (Displaying up to 5)\n", size(solutions,2));

num_sol = size(solutions,2);


% Ground truth

cam2_R = [ 0.9480719 , -0.01883705, -0.31749776 ;
           0.02188304,  0.99974235,  0.00602996 ;
           0.31730238, -0.01266465,  0.94823985 ];

cam2_t = [ 0.99828645,-0.05287068, 0.025077  ];

% cam2_aa: [-0.0095   -0.3230    0.0207];   % Converted from cam2_R

% Upside down R:
cam2_R2 = [ 0.95515593, -0.12487512, -0.26848343 ;
           -0.12268072, -0.99213115,  0.02500444 ;
           -0.26949321,  0.0090546 , -0.96295972 ];



for s = 1: min(num_sol, 5)
    R_err = StereoInterface.angles([0,0,1]*solutions(s).aa2mat()', [0,0,1] * cam2_R');
    for t = 1: min(numel(solutions(s).patches), 5)
    % Report R and t angular error
    t_err = StereoInterface.angles(solutions(s).patches(t).centre_xyz, cam2_t);
    fprintf("Solution %d:   aa = [%f %f %f], rpy: [%f %f %f], theta-phi = [%f %f], score: %d\n, R-error: %f (%f deg), t-error: %f (%f deg)\n", ...
            s, solutions(s).centre, RCube.R2rpy(solutions(s).aa2mat()), solutions(s).patches(t).centre, solutions(s).LB, ...
            R_err, rad2deg(R_err), t_err, rad2deg(t_err));
    
    R2_err = StereoInterface.angles([0,0,1]*solutions(s).aa2mat()', [0,0,1] * -cam2_R2');
    t2_err = StereoInterface.angles(solutions(s).patches(t).centre_xyz, -cam2_t);
    fprintf("-R2-error: %f (%f deg), -t-error: %f (%f deg)\n", R2_err, rad2deg(R2_err), t2_err, rad2deg(t2_err));
    % For each R, plot 1 solution matrix representing one possible t
    figure, imagesc(solutions(s).patches(t).edges_stop);
    end
end

%fprintf("Ground truth: aa = [%f %f %f], rpy = [%f %f %f], theta-phi = [%f %f]\n", ...
%        scene.cam2_aa, scene.cam2_rpy, scene.cam2_tp);

diary off;

% Save solutions to .mat file
save(['solutions\Rt_angel1_', expt_name, timestamp, '.mat'],'solutions');
