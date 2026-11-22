function [D,Sdown,info] = single_sphere_wall_stress_divergence( ...
    q,discretization,opt,step,checkHalfStep,stressletCheckResolution)
%SINGLE_SPHERE_WALL_STRESS_DIVERGENCE Height derivative of unit-force stresslet.
% q is one physical centre. opt.N_image is fixed for the entire stencil.
% D = d(S M e_z)/d q_z; lateral derivatives vanish by wall invariance.
% Sdown is the stresslet for a negative unit z-force, with zero torque.
% Uses 3 mobility solves, or 5 with the half-step check; no resistance solve.
% discretization has N, proxyRadius (fraction of a), collocationRatio.
% Rebuild all geometry at perturbed heights. Never difference across a change
% of source families/counts. q, step, stresses and derivatives are dimensional.

validateattributes(q,{'double'},{'real','finite','size',[1,3]});
validateattributes(step,{'double'},{'scalar','real','finite','positive'});
if q(3)-opt.wall_height-opt.radius <= step
    error('wall:StressStencil','The derivative stencil must remain above contact.');
end
if nargin < 6, stressletCheckResolution = []; end
offsets = [0,-step,step];
if checkHalfStep, offsets = [offsets,-step/2,step/2]; end
nSolves = numel(offsets);
stress = zeros(3,3,nSolves);
info = struct('offsets',offsets,'surfaceResidual',zeros(nSolves,1), ...
    'momentError',zeros(nSolves,1),'iterations',zeros(nSolves,1), ...
    'linearSystemResidual',zeros(nSolves,1),'converged',NaN(nSolves,1), ...
    'rank',zeros(nSolves,1),'stepChange',NaN,'stepChangeAbsolute',NaN, ...
    'stressletCorrection',NaN,'proxyCount',0,'sourceVectors',0, ...
    'collocationPoints',0,'imagePoints',0,'velocity',zeros(6,1));
for j = 1:nSolves
    centre = q+[0,0,offsets(j)];
    [rin,rout] = init_spheres(centre,discretization.proxyRadius, ...
        discretization.N,discretization.collocationRatio,opt.radius);
    load = [0;0;1;0;0;0];
    if j==1, load = -load; end
    [U,it,~,err,~,details] = solve_mobility_wall_images(centre,rin,rout,load,opt);
    stress(:,:,j) = wall_images_stresslet(details);
    basis = details.problem.body{1}.basis;
    counts = arrayfun(@(b) size(b.points,1),basis);
    if j==1
        baseCounts = counts; baseFamilies = [basis.family];
        info.proxyCount = details.N;
        info.sourceVectors = details.N_total;
        info.collocationPoints = details.M;
        info.velocity = U;
        for b = 1:numel(basis)
            if startsWith(basis(b).label,"wall ")
                info.imagePoints = counts(b); break
            end
        end
        if ~isempty(stressletCheckResolution)
            [~,check] = wall_images_stresslet(details,stressletCheckResolution);
            info.stressletCorrection = check.relativeCorrection;
        end
    elseif ~isequal(counts,baseCounts) || ~isequal([basis.family],baseFamilies)
        error('wall:StressStencil','Image basis changed across the derivative stencil.');
    end
    info.surfaceResidual(j) = err;
    info.momentError(j) = details.momentError;
    info.iterations(j) = it;
    info.linearSystemResidual(j) = details.linearSystemResidual;
    info.converged(j) = details.converged;
    info.rank(j) = sum(details.retainedRanks);
end
Sdown = stress(:,:,1);
D = (stress(:,:,3)-stress(:,:,2))/(2*step);
if checkHalfStep
    Dhalf = (stress(:,:,5)-stress(:,:,4))/step;
    info.stepChangeAbsolute = norm(Dhalf-D,'fro');
    info.stepChange = info.stepChangeAbsolute/max(norm(Dhalf,'fro'),realmin);
    D = Dhalf;
end
end
