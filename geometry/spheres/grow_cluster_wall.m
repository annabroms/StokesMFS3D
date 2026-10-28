function [q,geometry] = grow_cluster_wall(P,a,delta_body,delta_wall,seed,mode)
%GROW_CLUSTER_WALL Seeded equal-sphere cluster with explicit wall geometry mode.
% delta_body is an absolute surface gap, not a centre distance or gap ratio.
% Each inserted particle has a neighbour at that gap, and all gaps exceed it
% within tolerance. equalHeight changes only height during a wall-gap sweep.
% Returned wall is z=0; translate centres for a nonzero wall_height.
if nargin == 0, self_test(); q=[]; geometry=[]; return; end
if nargin < 6, mode = "equalHeight"; end
validateattributes(P,{'double'},{'scalar','integer','positive'});
validateattributes(a,{'double'},{'vector','positive','real','finite'});
assert((isscalar(a) || numel(a)==P) && all(a==a(1)), ...
    'wall:UnequalRadii','First version requires equal radii.');
a = a(1);
validateattributes(delta_body,{'double'},{'scalar','positive','finite'});
validateattributes(delta_wall,{'double'},{'scalar','positive','finite'});
validateattributes(seed,{'double'},{'scalar','integer','nonnegative','finite'});
mode = string(validatestring(mode,{'equalHeight','threeDimensional'}));
oldRng = rng; cleanup = onCleanup(@() rng(oldRng));
rng(seed,'twister');
tolerance = 1e-5*delta_body;
if mode == "equalHeight"
    lateral = a*grow_cluster(P,delta_body/a,2);
    q = [real(lateral),imag(lateral),(a+delta_wall)*ones(P,1)];
else
    q = a*grow_cluster(P,delta_body/a,3);
    q(:,3) = q(:,3)-min(q(:,3))+a+delta_wall;
end
nearest = inf(P,1);
for p = 1:P
    others = [1:p-1,p+1:P];
    if ~isempty(others), nearest(p)=min(vecnorm(q(others,:)-q(p,:),2,2)-2*a); end
end
bodyGap = NaN;
if P > 1
    bodyGap = min(nearest);
    assert(bodyGap > 0 && all(nearest>=delta_body-tolerance));
    assert(all(abs(nearest-delta_body)<=tolerance));
end
wallGap = min(q(:,3)-a);
assert(wallGap>0 && abs(wallGap-delta_wall)<=tolerance);
geometry = struct('P',P,'radius',a,'delta_body',delta_body, ...
    'delta_wall',delta_wall,'seed',seed,'mode',mode, ...
    'nearestBodyGaps',nearest,'minBodyGap',bodyGap,'minWallGap',wallGap, ...
    'tolerance',tolerance);
end

function self_test()
% Check reproducibility, singleton handling, both modes, and RNG restoration.
before = rng;
for mode = ["equalHeight","threeDimensional"]
    for P = [1,2,10]
        [q,g] = grow_cluster_wall(P,2,0.7,0.2,7,mode);
        assert(isequal(q,grow_cluster_wall(P,2,0.7,0.2,7,mode)));
        assert(abs(g.minWallGap-0.2)<g.tolerance);
        assert(size(q,1)==P);
        fprintf('grow_cluster_wall: PASS; %s, P=%d, minimum body gap %.6g, wall gap %.6g.\n', ...
            mode,P,g.minBodyGap,g.minWallGap);
        if P == 10
            plotOpt = struct('wall_height',0,'delta_body',0.7,'delta_wall',0.2, ...
                'visible',get(groot,'DefaultFigureVisible'));
            % Display the actual mode, including the distinct 3D geometry.
            wallOpt = wall_solver_options(struct('N',30,'backend',"blake")); wallOpt.radius = 2;
            [sources,~] = init_spheres(q,wallOpt.proxyRadius, ...
                wallOpt.N,wallOpt.collocationRatio,2);
            fig = plot_wall_geometry(q,2,sources,plotOpt);
            set(fig,'Name',char("grow_cluster_wall self-test: "+mode));
        end
    end
end
assert(isequal(before,rng));
fprintf('grow_cluster_wall: reproducibility and caller RNG restoration verified.\n');
end
