%DEMO_CLUSTER_RESISTANCE_WITH_IMAGES Ten close spheres in resistance motion.
% Date: 2026-10-16. Inspired by the paper's randomized P=20 cluster
% (Example 6 and Fig. 1), but uses P=10 for a readable direct-sum demo.
% Every sphere has a neighbour at gap about 1e-3. Random rigid velocities
% are prescribed; the solver returns the forces and torques needed to drive
% them. The scatter plot shows independent relative no-slip residuals.
%
% Set capMode to "fibonacci" for the paper's original two quasi-uniform
% caps per contact. "exponential" instead uses one cap per contact (first
% angle, first layer count and exponent; total point count of both caps)
% with points graded toward contact, and needs local areas.
% N=686 is the paper's coarse spherical design. N_image=[] chooses the
% paper's adaptive image count from each pair's gap (Eq. 39). This run can
% be costly because cross-sphere kernels are summed directly, without FMM.
% Reference: A. Broms, A. H. Barnett, and A.-K. Tornberg,
% "Accurate close interactions of Stokes spheres using lubrication-adapted
% image systems," Journal of Computational Physics 523 (2025), 113636.

clear;
close all;
fprintf(['Ten-sphere resistance demo: prescribed rigid velocities for a ', ...
    'close cluster; solve for forces and torques using gap-adaptive ', ...
    'image lines and collocation caps.\n']);
fprintf(['Reference: Broms, Barnett and Tornberg, "Accurate close ', ...
    'interactions of Stokes spheres using lubrication-adapted image ', ...
    'systems" (Journal of Computational Physics 523, 2025).\n']);
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));

P = 5;
gap = 1e-3;
capMode = "exponential"; % Change between exponential or "fibonacci" for the paper's caps.
rng(20261016);
[q,~] = grow_cluster(P,gap,3);
U = randn(6*P,1); % Per sphere: [translation; angular velocity].

if capMode == "exponential"
    rowWeightMode = "local";
else
    rowWeightMode = "paper";
end
N = 686;
N_image = []; % Adaptive, separately for each close pair.
opt = struct('N',N,'proxyRadius',0.63, ...
    'N_image',N_image,'imageFamilies',"SRD", ...
    'capMode',capMode,'capAngles',[pi/5,pi/60], ...
    'capUpsampling',[6,6],'capLayers',[18,18], ...
    'capExponent',[3,3],'rowWeightMode',rowWeightMode, ...
    'validationCount',180,'validationGapCount',301, ...
    'gmresTolerance',1e-6,'maxIterations',250);

[F,iters,lambdaNorm,maxResidual,details] = ...
    solve_resistance_with_images(q,U,opt);

fprintf('\nTen-sphere cluster: %d close pairs; cap mode %s\n', ...
    size(details.grid.pairs,1),capMode);
fprintf('N=%d; per-contact image counts from the gap rule:\n',details.grid.N);
fprintf('body pair     gap       N_image\n');
for pair = 1:size(details.grid.pairs,1)
    data = details.grid.pairs(pair,:);
    fprintf('%4d %4d   %.3e       %d\n',data);
end
fprintf('GMRES iterations %d, true weighted residual %.3e\n', ...
    iters,details.trueGmresResidual);
fprintf('Max source strength %.3e\n',lambdaNorm);
% Paper's measure: boundary-condition error over max prescribed speed.
pointError = cell(P,1);
maxSpeed = 0;
for k = 1:P
    points = details.validationPoints{k};
    desired = reshape(getKmat(points,q(k,:))*U(6*(k-1)+(1:6)),3,[]).';
    pointError{k} = vecnorm(details.validationVelocity{k}-desired,2,2);
    maxSpeed = max(maxSpeed,max(vecnorm(desired,2,2)));
end
pointError = cellfun(@(e) e/maxSpeed,pointError,'UniformOutput',false);
bcError = cellfun(@max,pointError);
maxBC = max(bcError);
fprintf('Max boundary-condition error relative to max velocity %.3e\n',maxBC);
fprintf('Body       Fx       Fy       Fz       Tx       Ty       Tz    max BC\n');
for k = 1:P
    ft = F(6*(k-1)+(1:6));
    fprintf('%4d %8.2g %8.2g %8.2g %8.2g %8.2g %8.2g   %.2e\n', ...
        k,ft,bcError(k));
end

fig = figure('Name','Ten-sphere cluster: surface residual');
ax = axes(fig);
hold(ax,'on');
for k = 1:P
    points = details.validationPoints{k};
    scatter3(ax,points(:,1),points(:,2),points(:,3), ...
        8,log10(max(pointError{k},eps)),'filled');
end
axis(ax,'equal');
grid(ax,'on');
view(ax,3);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
colormap(ax,turbo);
cb = colorbar(ax);
cb.Label.String = 'log_{10} BC error / max velocity';
title(ax,sprintf('P=%d, gap=%.0e, %s caps; max BC error %.2e', ...
    P,gap,capMode,maxBC));
