function [C,W,KR] = ring_constant_maps(ringId)
%RING_CONSTANT_MAPS Ring-to-node copy map and force maps of note Section 10.
%
%   [C,W,KR] = ring_constant_maps(ringId)
%
% Input:  ringId N-by-1 ring index of each source node (1..R).
% Outputs (all Cartesian, ring/node-interleaved [x1;y1;z1;x2;...]):
%   C  3N-by-3R sparse copy map, C = L (x) I3, lambda = C*a, eq. (64).
%      Every nonzero block is I3; nothing rotates with azimuth.
%   W  3R-by-3, W = [n_1 I3; ...; n_R I3]; W'*a is the total force of
%      per-node ring strengths a (3R-by-1), eq. (67)-(68).
%   KR 3R-by-3, KR = ones(R,1) (x) I3, so W'*KR = N*I3.
if nargin == 0
    self_test();
    return
end
ringId = ringId(:);
N = numel(ringId);
R = max(ringId);
nr = accumarray(ringId,1);
L = sparse((1:N).',ringId,1,N,R);
C = kron(L,speye(3));
W = kron(nr,eye(3));
KR = kron(ones(R,1),eye(3));
end

function self_test()
%SELF_TEST Check eq. (64), (68) on a small ring layout.
ringId = [1;1;1;2;2;2;2;3;3;3];
N = numel(ringId); R = 3;
[C,W,KR] = ring_constant_maps(ringId);
assert(isequal(size(C),[3*N,3*R]) && isequal(size(W),[3*R,3]));
a = reshape(1:9,3,3);                   % columns are ring vectors a_r
lambda = reshape(C*a(:),3,[]);
assert(isequal(lambda(:,4),a(:,2)) && isequal(lambda(:,10),a(:,3)));
assert(norm(C.'*kron(ones(N,1),eye(3))-W,'fro')==0);   % K_N'*C = W'
assert(norm(W.'*KR-N*eye(3),'fro')==0);
assert(norm(W.'*a(:)-sum(lambda,2))<1e-13);
fprintf('ring_constant_maps self-test passed.\n');
end
