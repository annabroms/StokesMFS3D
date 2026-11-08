function [alpha,a,U,info] = ring_complete_solve(A,V,W,KR,a0,N,beta,rowScale,colScale,tol)
%RING_COMPLETE_SOLVE Force-completed ring-coefficient solve, note eqs. (74)-(79).
%
%   [alpha,a,U,info] = ring_complete_solve(A,V,W,KR,a0,N,beta,rowScale,colScale,tol)
%
% Solves  (A*Pi_R + beta*V*W') alpha = -A*a0,  Pi_R = I - KR*W'/N,
% using the rank-three form  A + (beta*V - A*KR/N)*W'  (eq. 79); the dense
% projector is never formed. Then a = Pi_R*alpha + a0 and U = -beta*W'*alpha.
% A*a - V*U = 0 are the boundary rows, so W'*a = Fg holds exactly.
%
% Inputs:  A (3M-by-3R) boundary rows applied to per-node ring vectors
%            a = [a_1;...;a_R] (3R-by-1, ring-major [x;y;z]);
%          V (3M-by-3) translation coupling; W,KR (3R-by-3) from
%            ring_constant_maps; a0 (3R-by-1) completion KR*Fg/N; N sources;
%          beta velocity/force conversion scale (not the mobility);
%          rowScale (3M-by-1) multiplies rows; colScale (scalar or 3R-by-1)
%            gives alpha = colScale.*alpha_hat; tol relative SVD cutoff.
% Outputs: alpha (3R-by-1) auxiliary coefficients (not physical strengths),
%          a (3R-by-1) physical per-node ring vectors, U (3-by-1) translation,
%          info: size, sigmaMax, threshold (= tol*sigmaMax), rank.
if nargin == 0
    self_test();
    return
end
Msys = A+(beta*V-(A*KR)/N)*W.';
Msys = rowScale(:).*Msys.*(colScale(:).');
rhs = -rowScale(:).*(A*a0);
info.size = size(Msys);
fprintf('Ring system matrix %d x %d before SVD.\n',info.size);
[Y,Uleft,sigma] = getPseudoFactors(Msys,tol,0);
alpha = colScale(:).*(Y*(Uleft.'*rhs));
t = W.'*alpha;                      % force-weighted ring sum, eq. (77)
a = alpha-KR*t/N+a0;
U = -beta*t;
info.sigmaMax = sigma(1);
info.threshold = tol*info.sigmaMax;
info.rank = numel(sigma);
fprintf('Truncation threshold %.3e (relative %.1e); numerical rank %d of %d.\n', ...
    info.threshold,tol,info.rank,min(info.size));
end

function self_test()
%SELF_TEST Recover a known (a,U); check eq. (79) and force exactness under truncation.
rng(3);
R = 5; N = 23; nr = [3;5;6;5;4]; M = 40;
[~,W,KR] = ring_constant_maps(repelem((1:R).',nr));
Fg = [0;0;-2]; beta = 0.7; a0 = KR*Fg/N;
Pi = eye(3*R)-KR*W.'/N;
assert(norm(W.'*Pi,'fro')<1e-13);
aTrue = a0+Pi*randn(3*R,1); UTrue = randn(3,1);
A = randn(3*M,3*R); V = randn(3*M,3);
% Adjust V by a rank-one term so the known pair satisfies A*a-V*U = 0.
V = V+(A*aTrue-V*UTrue)*(UTrue.'/(UTrue.'*UTrue));
assert(norm(A*aTrue-V*UTrue)<1e-10);
alphaTrue = (aTrue-a0)-KR*UTrue/(beta*N);       % eq. (76)
assert(norm(Pi*alphaTrue+a0-aTrue)<1e-12 && norm(-beta*W.'*alphaTrue-UTrue)<1e-12);
assert(norm((A*Pi+beta*V*W.')-(A+(beta*V-A*KR/N)*W.'),'fro')<1e-10);   % eq. (79)
[~,a,U,info] = ring_complete_solve(A,V,W,KR,a0,N,beta,ones(3*M,1),1,1e-12);
assert(norm(a-aTrue)<1e-8*norm(aTrue) && norm(U-UTrue)<1e-8*norm(UTrue));
assert(info.rank<=3*R);
[~,a2,~,info2] = ring_complete_solve(A,V,W,KR,a0,N,beta,ones(3*M,1),1,0.3);
assert(info2.rank<info.rank && norm(W.'*a2-Fg)<1e-12);   % force exact after truncation
fprintf('ring_complete_solve self-test passed.\n');
end
