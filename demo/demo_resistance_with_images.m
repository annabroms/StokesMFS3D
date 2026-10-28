%DEMO_RESISTANCE_WITH_IMAGES Small two-sphere resistance solve.
% Run setpath first, then this demo.
% Prescribed opposite translations cause a lubrication interaction.
% Uses the paper's N=686 coarse sources and adaptive N_image per contact.
% Reference: A. Broms, A. H. Barnett, and A.-K. Tornberg,
% "Accurate close interactions of Stokes spheres using lubrication-adapted
% image systems," Journal of Computational Physics 523 (2025), 113636.

useSVD = true; % true: dense global SVD solve; false: one-body preconditioned GMRES.

fprintf(['Two-sphere resistance demo: prescribe opposite translations ', ...
    'across a narrow gap, then solve for forces and torques using ', ...
    'gap-adaptive image lines.\n']);
fprintf(['Reference: Broms, Barnett and Tornberg, "Accurate close ', ...
    'interactions of Stokes spheres using lubrication-adapted image ', ...
    'systems" (Journal of Computational Physics 523, 2025).\n']);
q = [0,0,0;2.05,0,0];
U = [1;0;0;0;0;0;-1;0;0;0;0;0];
N = 686;
N_image = []; % Adaptive from the sphere gap by paper Eq. (39).
capMode = "exponential"; % Change to "fibonacci" for the paper's caps.
capPointCounts = [];
if capMode == "exponential"
    % Total cap nodes per contact, overriding the image-based count
    % (6*3*N_image per cap, only 108 at this gap), which is too sparse for
    % 18 exponentially graded rings.
    capPointCounts = 1078;
end
opt = struct('N',N,'N_image',[], ...
    'capMode',capMode,'capPointCounts',capPointCounts, ...
    'capLayers',[18,18],'capExponent',[3,3],'capMinAzimuth',8, ...
    'validationCount',1000, ...
    'validationGapCount',401);
if useSVD
    opt.solverMode = "denseSVD";
end
[Fvec,iters,lambdaNorm,maxResidual,details] = ...
    solve_resistance_images(q,U,opt);
if useSVD
    fprintf('Dense SVD solve, rank %d, linear system residual: %.3e\n', ...
        details.retainedRanks,details.linearSystemResidual);
else
    fprintf('Iterations: %d, true weighted residual: %.3e\n', ...
        iters,details.trueGmresResidual);
end
fprintf('N=%d, N_image=%d per close contact\n', ...
    details.grid.N,details.grid.pairs(1,4));
if ~isempty(capPointCounts)
    fprintf('Requested exponential cap nodes per contact: %d\n', ...
        capPointCounts);
end
cap = details.grid.body(1).capPoints;
dir = (q(2,:)-q(1,:))/norm(q(2,:)-q(1,:));
theta = acos(max(-1,min(1,(cap-q(1,:))*dir.')));
fprintf('%s caps, %s row areas: %d cap nodes on sphere 1\n', ...
    details.options.capMode,details.options.rowWeightMode,size(cap,1));
fprintf('Nearest cap node: %.3e rad; %d within quarter of inner cap\n', ...
    min(theta),nnz(theta <= details.options.capAngles(2)/4));
fprintf('Max independent boundary residual: %.3e\n',maxResidual);
fprintf('Max source strength: %.3e\n',lambdaNorm);
disp('Per-sphere [force; torque] (existing source-strength sign):');
disp(reshape(Fvec,6,[]).');

fig = figure('Name','Image resistance demo: sphere 1');
fig.Position(3:4) = [1400,460];
layout = tiledlayout(fig,1,3,'TileSpacing','compact');
ax = nexttile(layout);
points = details.validationPoints{1};
scatter3(ax,points(:,1),points(:,2),points(:,3),8, ...
    log10(max(details.validationRelative{1},eps)),'filled');
axis(ax,'equal'); view(ax,3); colorbar(ax);
xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
title(ax,'log_{10} independent relative surface residual');

% The spheres lie on the x axis, so these are actual spherical-cap
% coordinates, with x showing the surface curvature toward contact.
for panel = 1:2
    a = details.options.capAngles(panel);
    onCap = theta <= a+10*eps;
    [az,pol] = meshgrid(linspace(0,2*pi,60),linspace(0,a,24));
    x = q(1,1)+cos(pol);
    y = q(1,2)+sin(pol).*cos(az);
    z = q(1,3)+sin(pol).*sin(az);
    ax = nexttile(layout);
    surf(ax,x,y,z,'FaceColor',[0.75,0.8,0.85], ...
        'FaceAlpha',0.12,'EdgeColor','none');
    hold(ax,'on');
    scatter3(ax,cap(onCap,1),cap(onCap,2),cap(onCap,3), ...
        10,[0.1,0.35,0.85],'filled');
    edge = linspace(0,2*pi,200);
    plot3(ax,q(1,1)+cos(a)*ones(size(edge)), ...
        q(1,2)+sin(a)*cos(edge),q(1,3)+sin(a)*sin(edge), ...
        'r-','LineWidth',1.2);
    scatter3(ax,q(1,1)+1,q(1,2),q(1,3),35,'k','filled');
    pad = 0.15*(1-cos(a));
    lim = 1.1*sin(a);
    xlim(ax,q(1,1)+[cos(a)-pad,1+pad]);
    ylim(ax,q(1,2)+[-lim,lim]);
    zlim(ax,q(1,3)+[-lim,lim]);
    axis(ax,'equal'); grid(ax,'on');
    view(ax,[1,0.5,0.4]);
    xlabel(ax,'x'); ylabel(ax,'y'); zlabel(ax,'z');
    if panel == 1
        title(ax,sprintf('%s cap points on sphere 1',capMode));
    else
        title(ax,'Contact close-up on the sphere');
    end
end
