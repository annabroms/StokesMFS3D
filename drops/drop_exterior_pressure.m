function pressure = drop_exterior_pressure( ...
    targets, innerSources, exteriorStrengths, targetBatchSize)
%DROP_EXTERIOR_PRESSURE Sum exterior pressure from all drop Stokeslets.
%
% INPUTS
%   targets           Q-by-3 exterior evaluation points.
%   innerSources      One cell of N-by-3 inner proxy points per drop.
%   exteriorStrengths One cell of interleaved 3N-by-1 force-free lambda
%                     coefficients per drop.
%   targetBatchSize   Maximum target count per direct kernel call. This
%                     limits temporary memory; it does not alter accuracy.
%
% OUTPUT
%   pressure          Q-by-1 exterior pressure, with p -> 0 at infinity.
%
% The drop model has p_0(x)=sum_i sum_j P(x,y_inner_ij)*lambda_ij.
% getPressure supplies P and its 1/(4*pi) normalization. This routine only
% batches the targets and adds contributions from different drops.

if nargin == 0
    self_test();
    return
end

pointCount = size(targets,1);
pressure = zeros(pointCount,1);
for first = 1:targetBatchSize:pointCount
    selected = first:min(first+targetBatchSize-1,pointCount);
    for body = 1:numel(innerSources)
        pressure(selected) = pressure(selected) + getPressure( ...
            exteriorStrengths{body},innerSources{body},targets(selected,:));
    end
end
end

function self_test()
% Verify two-drop superposition and identical answers for different batches.
innerSources = {[-2 0 0;-1.5 0 0],[1.5 0 0;2 0 0]};
exteriorStrengths = {[1;0;0;-1;0;0],[-1;0;0;1;0;0]};
targets = [0 0 0;0 0.5 0;0 0 1;0.2 0.3 0.4];
reference = getPressure( ...
    [exteriorStrengths{1};exteriorStrengths{2}], ...
    [innerSources{1};innerSources{2}],targets);
oneAtATime = drop_exterior_pressure( ...
    targets,innerSources,exteriorStrengths,1);
allAtOnce = drop_exterior_pressure( ...
    targets,innerSources,exteriorStrengths,size(targets,1));
if ~isequal(size(allAtOnce),[size(targets,1),1]) || ...
        norm(reference-oneAtATime,inf)>1e-14 || ...
        norm(reference-allAtOnce,inf)>1e-14
    error('drop_exterior_pressure self-test failed.');
end
fprintf('drop_exterior_pressure self-test passed.\n');
end
