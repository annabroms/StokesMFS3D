function pressure = getPressure(strengths, sources, targets)
%GETPRESSURE Evaluate pressure from three-dimensional Stokeslet sources.
%
%   pressure = getPressure(strengths, sources, targets)
%
% INPUTS
%   strengths  3N-by-1 source strengths, interleaved by source node as
%              [fx1;fy1;fz1;fx2;fy2;fz2;...].
%   sources    N-by-3 Stokeslet locations.
%   targets    M-by-3 pressure evaluation locations.
%
% OUTPUT
%   pressure   M-by-1 pressure, with zero pressure at infinity.
%
% For r = target-source, one strength f contributes (r dot f)/(4*pi*|r|^3).
% The pressure kernel has no viscosity factor. This direct implementation
% matches the Stokeslet normalization in stokes_SLP_mat/getStokesletFlow.
% A source must not coincide with a target; batch targets in the caller
% when the M-by-N pair arrays would be too large.

if nargin == 0
    self_test();
    return
end

force = reshape(strengths,3,[]).';
dx = targets(:,1)-sources(:,1).';
dy = targets(:,2)-sources(:,2).';
dz = targets(:,3)-sources(:,3).';
distanceSquared = dx.^2+dy.^2+dz.^2;
if any(distanceSquared==0,'all')
    error('getPressure:CoincidentPoints', ...
        'A pressure target coincides with a Stokeslet source.');
end
radialForce = dx.*force(:,1).' + dy.*force(:,2).' + ...
    dz.*force(:,3).';
pressure = sum(radialForce./(4*pi*distanceSquared.^(3/2)),2);
end

function self_test()
% Check the signed reference value, kernel normalization, and superposition.
targets = [1 0 0;-1 0 0;0 1 0];
actual = getPressure([1;0;0],[0 0 0],targets);
expected = [1;-1;0]/(4*pi);
if norm(actual-expected,inf)>1e-14
    error('getPressure self-test: single Stokeslet reference failed.');
end

sources = [-1 0 0;1 0 0];
strengths = [1;0;0;-1;0;0];
targets = [0 0 0;0 1 0;2 0 0];
combined = getPressure(strengths,sources,targets);
separate = getPressure(strengths(1:3),sources(1,:),targets) + ...
    getPressure(strengths(4:6),sources(2,:),targets);
if norm(combined-separate,inf)>1e-14 || ...
        abs(combined(1)-1/(2*pi))>1e-14
    error('getPressure self-test: superposition or prefactor failed.');
end
fprintf('getPressure self-test passed.\n');
end
