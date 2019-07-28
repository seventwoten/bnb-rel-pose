% Sample scene set-up:
% Create an array s30, containing 10 stereo scenes with 30 points each.
% (Replace 30 with any number of points)

n_points = 30;
n_scenes = 10;

s30 = Scene.empty;
cmap = parula(n_points);

for i = 1:n_scenes
    s30(i) = Scene([], n_points);
    s30(i) = s30(i).setCam2();
    s30(i) = s30(i).setView1();
    s30(i) = s30(i).setView2();
    
    % Show pairs of views
    figure;
    subplot(1,2,1), scatter(s30(i).view1(:,1), s30(i).view1(:,2), [], cmap), hold on, xlim([-1,1]), ylim([-1,1]);
    subplot(1,2,2), scatter(s30(i).view2(:,1), s30(i).view2(:,2), [], cmap), hold on, xlim([-1,1]), ylim([-1,1]);    
end

save('scenes\s30.mat','s30');