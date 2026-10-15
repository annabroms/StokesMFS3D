function [maximumSeparation, rmsSeparation, lengthDifference] = ...
    compare_streamline_paths(coarsePaths, finePaths, seeds, sampleCount)
%COMPARE_STREAMLINE_PATHS Compare matched bidirectional x-z streamlines.
%
% INPUTS
%   coarsePaths, finePaths  Cell arrays from get_streamline, with the same
%                           seed order and the same tracing domain.
%   seeds                    K-by-2 starting positions [x,z].
%   sampleCount              Number of equally spaced arc-length samples
%                           in EACH direction from every seed.
%
% OUTPUTS (one K-by-1 value per seed)
%   maximumSeparation        Largest point-to-point distance over both
%                           directions at equal distance from the seed.
%   rmsSeparation            Root-mean-square of those distances.
%   lengthDifference         Largest absolute difference between traced
%                           forward and backward arc lengths.
%
% Paths are sampled only as far as both computations traced. This avoids
% comparing an endpoint of one path with unrelated positions beyond the
% endpoint of the other. lengthDifference reports the omitted tail size.
% If one trace consists only of its seed, the common-path separation is
% zero and lengthDifference records how far the other trace travelled.

if nargin == 0
    self_test();
    return
end

pathCount = numel(coarsePaths);
maximumSeparation = zeros(pathCount,1);
rmsSeparation = zeros(pathCount,1);
lengthDifference = zeros(pathCount,1);
for pathIndex = 1:pathCount
    coarse = coarsePaths{pathIndex};
    fine = finePaths{pathIndex};
    seed = seeds(pathIndex,:);
    [~,coarseSeedIndex] = min(vecnorm(coarse-seed,2,2));
    [~,fineSeedIndex] = min(vecnorm(fine-seed,2,2));
    coarseDirections = {flipud(coarse(1:coarseSeedIndex,:)), ...
        coarse(coarseSeedIndex:end,:)};
    fineDirections = {flipud(fine(1:fineSeedIndex,:)), ...
        fine(fineSeedIndex:end,:)};
    allDistances = [];
    for direction = 1:2
        coarseTrace = coarseDirections{direction};
        fineTrace = fineDirections{direction};
        % Specify the row dimension: diff on a single [x,z] row would
        % otherwise subtract its two coordinates and invent a path step.
        coarseArc = [0;cumsum(vecnorm(diff(coarseTrace,1,1),2,2))];
        fineArc = [0;cumsum(vecnorm(diff(fineTrace,1,1),2,2))];
        lengthDifference(pathIndex) = max(lengthDifference(pathIndex), ...
            abs(coarseArc(end)-fineArc(end)));
        commonLength = min(coarseArc(end),fineArc(end));
        if commonLength == 0
            % A stagnant or immediately stopped path has only its seed.
            distances = norm(coarseTrace(1,:)-fineTrace(1,:));
        else
            arcSamples = linspace(0,commonLength,sampleCount).';
            [coarseArc,coarseUnique] = unique(coarseArc,'stable');
            [fineArc,fineUnique] = unique(fineArc,'stable');
            coarseSamples = interp1(coarseArc,coarseTrace(coarseUnique,:), ...
                arcSamples,'linear');
            fineSamples = interp1(fineArc,fineTrace(fineUnique,:), ...
                arcSamples,'linear');
            distances = vecnorm(coarseSamples-fineSamples,2,2);
        end
        allDistances = [allDistances;distances]; %#ok<AGROW>
    end
    maximumSeparation(pathIndex) = max(allDistances);
    rmsSeparation(pathIndex) = sqrt(mean(allDistances.^2));
end
end

function self_test()
% Exercise different ode-like node spacings, a shifted line, and a
% single-point path caused by immediate stopping near stagnation.
seeds = [0 0;0 0;0 0];
coarse = {[-1 0;0 0;1 0],[-1 0;0 0;1 0],[0 0]};
fine = {[-1 0;-0.5 0;0 0;0.25 0;1 0], ...
    [-1 0.1;-0.5 0.1;0 0;0.5 0.1;1 0.1], ...
    [0 0;0.5 0;1 0]};
[maximum,rmsValue,lengthGap] = compare_streamline_paths( ...
    coarse,fine,seeds,21);
if maximum(1)>1e-12 || rmsValue(1)>1e-12 || ...
        lengthGap(1)>1e-12 || maximum(2)<0.05 || ...
        maximum(3)>1e-12 || rmsValue(3)>1e-12 || ...
        abs(lengthGap(3)-1)>1e-12
    error('compare_streamline_paths self-test failed.');
end
fprintf('compare_streamline_paths self-test passed.\n');
end
