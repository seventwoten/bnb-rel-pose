function vis_sol(solutions, scene, fignum)
    if ~exist('fignum', 'var') || isempty(fignum)
        fignum = 6;
    end

    figure(fignum);
    figure(fignum+1), view(3), xlim([-1.5 1.5]), ylim([-1.5 1.5]), zlim([-1.5 1.5]), ...
        xlabel('x'), ylabel('y'), zlabel('z'), hold on;
    cmap = parula(numel(solutions));

    r = zeros(numel(solutions), 3);
    c = cell(1,numel(solutions));
    for i = 1:numel(solutions)
        r(i,:) = solutions(i).centre;
        c{i} = zeros(numel(solutions(i).patches), 3);
        for t = 1:numel(solutions(i).patches)
            c{i}(t, :) = solutions(i).patches(t).centre_xyz;
        end
    end

    figure(fignum), scatter3(r(:,1),r(:,2),r(:,3),20,cmap,'filled'), ...
        xlim([-pi pi]), ylim([-pi pi]), zlim([-pi pi]),...
        xlabel('x'), ylabel('y'), zlabel('z'), hold on; 
    
    for i = 1:numel(solutions)
        figure(fignum+1), scatter3(c{i}(:,1),c{i}(:,2),c{i}(:,3),20,cmap(i,:),'filled'); 
        hold on;
    end
    
    % Mark ground truth with crosses
    if exist('scene', 'var') && ~isempty(scene)
        figure(fignum),   scatter3(scene.cam2_aa(1),scene.cam2_aa(2),scene.cam2_aa(3),40,'rx');
        t_xyz = tPatch.spherical2Cartesian(1, scene.cam2_tp(1), scene.cam2_tp(2)); 
        figure(fignum+1), scatter3(t_xyz(1),t_xyz(2),t_xyz(3),40,'rx');
    end
end