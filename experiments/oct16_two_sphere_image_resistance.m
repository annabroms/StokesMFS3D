%OCT16_TWO_SPHERE_IMAGE_RESISTANCE Reproduce the paper's Fig. 10 cases.
% Run this script for selected cases. Set runFullSweep below to true for
% the convergence sweep. Results remain in the workspace.
% Date: 2026-10-16. Execute the full sweep deliberately.
% Reference: A. Broms, A. H. Barnett, and A.-K. Tornberg,
% "Accurate close interactions of Stokes spheres using lubrication-adapted
% image systems," Journal of Computational Physics 523 (2025), 113636.
%
% Two unit spheres squeeze, shear, or counter-rotate at fixed gaps. All
% enhanced runs use N=686 spherical-design proxy nodes at Rp=0.63, two
% paper Fibonacci caps per contact, S/R/D image strengths, paper row
% weights and relative SVD cutoff 5e-12. The solverMode setting below
% selects one-body-preconditioned GMRES or a dense weighted SVD. Squeezing
% force is compared with Eq. (32), the Brenner series. Independent residual
% nodes include a 5000-request design plus 2001 near-gap cap points. Run
% only the final Plot results section to redraw saved results.

%% Run resistance experiments
clear;
close all;
runFullSweep = false; % Set true for the full image-count sweep.
solverMode = "oneBody"; % Choose "oneBody" or "denseSVD".
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root));
gaps = [1e-2,1e-3];
N_image_values = [0,10,20,30];
if runFullSweep
    gaps = logspace(-3,-1,8);
    N_image_values = 0:2:40;
end
resultTemplate = struct('motion',"",'gap',0,'N_image',0, ...
    'solverMode',solverMode, ...
    'force',zeros(12,1),'brennerForce',NaN,'forceError',NaN, ...
    'maxResidual',0,'iterations',0,'trueGmresResidual',0, ...
    'linearSystemResidual',0,'coefficientNorm',0,'retainedRanks',[]);
results = repmat(resultTemplate,3*numel(gaps)*numel(N_image_values),1);
resultIndex = 0;
motions = ["squeeze","shear","counter-rotation"];
for motionIndex = 1:3
    for gap = gaps
        q = [0,0,0;2+gap,0,0];
        U = zeros(12,1);
        switch motionIndex
            case 1
                U([1,7]) = [1,-1];
            case 2
                U([2,8]) = [1,-1];
            case 3
                U([6,12]) = [1,-1];
        end
        reference = NaN;
        if motionIndex == 1
            % Brenner's Eq. (32) for unit-speed squeezing, viscosity one.
            alpha = acosh(1+gap/2);
            series = 0;
            for n = 1:2000
                numerator = 4*cosh((n+0.5)*alpha)^2+ ...
                    (2*n+1)^2*sinh(alpha)^2;
                denominator = 2*sinh((2*n+1)*alpha)- ...
                    (2*n+1)*sinh(2*alpha);
                term = n*(n+1)/((2*n-1)*(2*n+3))* ...
                    (numerator/denominator-1);
                series = series+term;
                if n > 10 && abs(term) < 1e-12*abs(series)
                    break
                end
            end
            reference = 8*pi*sinh(alpha)*series;
        end
        for index = 1:numel(N_image_values)
            N_image = N_image_values(index);
            opt = struct('N',686,'proxyRadius',0.63, ...
                'N_image',N_image,'imageFamilies',"SRD", ...
                'capMode',"fibonacci",'capAngles',[pi/5,pi/60], ...
                'capUpsampling',[6,6],'rowWeightMode',"paper", ...
                'svdCutoff',5e-12,'gmresTolerance',1e-6, ...
                'solverMode',solverMode, ...
                'validationCount',5000,'validationGapCount',2001);
            [F,it,lambdaNorm,err,detail] = ...
                solve_resistance_with_images(q,U,opt);
            forceError = NaN;
            if motionIndex == 1
                forceError = abs(norm(F(1:3))-reference)/reference;
            end
            resultIndex = resultIndex+1;
            results(resultIndex) = struct('motion',motions(motionIndex), ...
                'gap',gap,'N_image',N_image, ...
                'solverMode',solverMode,'force',F, ...
                'brennerForce',reference,'forceError',forceError, ...
                'maxResidual',err,'iterations',it, ...
                'trueGmresResidual',detail.trueGmresResidual, ...
                'linearSystemResidual',detail.linearSystemResidual, ...
                'coefficientNorm',lambdaNorm, ...
                'retainedRanks',detail.retainedRanks);
            fprintf(['%s, %s gap=%g N_image=%d: BC %.3e, force error %.3e, ', ...
                'linear residual %.3e, it %d, |lambda| %.3e\n'], ...
                solverMode,motions(motionIndex),gap,N_image, ...
                err,forceError,detail.linearSystemResidual,it,lambdaNorm);
        end
    end
end

%% Plot results (run this section alone to replot existing results)
assert(exist('results','var') == 1 && ~isempty(results), ...
    'Run the experiment section or load a saved results variable first.');
figure('Name','2026-10-16 two-sphere image resistance');
layout = tiledlayout(1,3);
motionNames = ["squeeze","shear","counter-rotation"];
for motionIndex = 1:numel(motionNames)
    ax = nexttile(layout,motionIndex);
    hold(ax,'on'); grid(ax,'on');
    for gapValue = unique([results.gap])
        selected = results([results.motion] == motionNames(motionIndex) & ...
            [results.gap] == gapValue);
        if isempty(selected)
            continue
        end
        [counts,order] = sort([selected.N_image]);
        residuals = [selected.maxResidual];
        semilogy(ax,counts,residuals(order),'o-','DisplayName', ...
            sprintf('gap %.1e',gapValue));
    end
    yline(ax,1e-3,'k--','10^{-3} target');
    xlabel(ax,'N_image per close contact');
    ylabel(ax,'Maximum independent relative residual');
    title(ax,motionNames(motionIndex));
    legend(ax,'Location','best');
end
