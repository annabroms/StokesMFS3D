function [rinUp,routUp,transVec] = updateNodes(rin,rout,P,N,M,vel)
%UPDATENODES Translate proxy and collocation nodes by rigid-body velocities.
%
%   [rinUp,routUp,transVec] = UPDATENODES(rin,rout,P,N,M,vel)
%
%   Applies the translational components of vel to all nodes for each body.
%
%   rin, rout  - (N*P) x 3 and (M*P) x 3 node arrays (proxy, collocation)
%   P, N, M    - number of bodies, nodes per body (proxy, collocation)
%   vel        - 6P x 1 rigid body velocity (translation + rotation)
%
% TODO: account for rotations of the bodies

transVec = zeros(P,3);
rinUp = zeros(N*P,3);
routUp = zeros(M*P,3);

for k = 1:P
    velTrans = vel(6*(k-1)+1:6*(k-1)+3)';
    transVec(k,:) = velTrans;
    rinUp(N*(k-1)+1:k*N,:) = rin(N*(k-1)+1:k*N,:)+velTrans;
    routUp(M*(k-1)+1:k*M,:) = rout(M*(k-1)+1:k*M,:)+velTrans;
end
end
