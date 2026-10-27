function [rin,rout,opt] = init_spheres_wall(q,a,opt)
%INIT_SPHERES_WALL Standard spherical proxy sources, with no enhancement.
if nargin == 0, self_test(); rin = []; rout = []; opt = []; return; end
if nargin < 3, opt = struct(); end
opt.radius = a;
opt = wall_solver_options(opt);
assert(isscalar(a) || numel(a) == size(q,1));
a = a(1);
validateattributes(q,{'double'},{'real','finite','2d','ncols',3,'nonempty'});
assert(all(q(:,3)-a > opt.wall_height),'wall:WallIntersection', ...
    'Every particle must lie strictly above the wall.');
gridOpt = struct('des_n',opt.N,'a_glob',opt.collocationRatio);
[source,target] = getDesignGrid(opt.proxyRadius,gridOpt);
rin = zeros(size(q,1)*size(source,1),3);
rout = zeros(size(q,1)*size(target,1),3);
for p = 1:size(q,1)
    rin((p-1)*size(source,1)+(1:size(source,1)),:) = a*source+q(p,:);
    rout((p-1)*size(target,1)+(1:size(target,1)),:) = a*target+q(p,:);
end
end

function self_test()
% Check radius scaling and standard proxy placement without solving.
[s,t,o] = init_spheres_wall([0,0,3],2,struct('N',30,'backend',"blake"));
assert(max(abs(vecnorm(s-[0,0,3],2,2)-2*o.proxyRadius)) < 1e-12);
assert(max(abs(vecnorm(t-[0,0,3],2,2)-2)) < 1e-12);
fprintf('init_spheres_wall: PASS; a=2, %d sources at radius %.3g, %d collocation points.\n', ...
    size(s,1),2*o.proxyRadius,size(t,1));
figure('Name','Wall sphere grids: self-test');
scatter3(t(:,1),t(:,2),t(:,3),30,[0.2,0.5,0.8],'filled'); hold on;
scatter3(s(:,1),s(:,2),s(:,3),30,[0.8,0.3,0.2],'filled');
axis equal; grid on; xlabel('x'); ylabel('y'); zlabel('z'); view(35,25);
legend('Surface collocation','Interior proxy sources','Location','best');
title(sprintf('a=2, Rp/a=%.2f: standard proxy and collocation grids',o.proxyRadius));
end
