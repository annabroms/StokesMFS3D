function [RFD, U1, U2] = compute_RFD(q,rvec_in,rvec_out,RandVel,opt,precond)
%COMPUTE_RFD Random finite difference (RFD) drift term.
%
%   [RFD, U1, U2] = COMPUTE_RFD(q,rvec_in,rvec_out,RandVel,opt,precond)
%
%   If iterate_RFD is true, computes RFD using iterative solves via
%   solve_brownian_mobility. Otherwise, computes RFD using dense solves (for debugging).
%
%   Inputs:
%     q          - P x 3 particle centers
%     rvec_in    - (N*P) x 3 proxy nodes
%     rvec_out   - (M*P) x 3 collocation nodes
%     RandVel    - 6P x 1 random rigid body velocity
%     opt        - options struct:
%                  * opt.P, opt.N, opt.M (optional; inferred if missing)
%                  * opt.delta (required)
%                  * opt.iterate_RFD (required)
%                  * opt.gmres_tol, opt.fmm, opt.dt, etc.
%     precond    - struct with fields: Y, UU, LL, Kin, Tblock, Ktot, Btot
%
%   Outputs:
%     RFD - 6P x 1 finite difference drift term
%     U1, U2 - rigid body velocities for down/up configurations

if nargin < 5
    precond = [];
end

if ~isfield(opt,'delta')
    error('compute_RFD: opt.delta required');
end
if ~isfield(opt,'iterate_RFD')
    error('compute_RFD: opt.iterate_RFD required');
end

if isfield(opt,'P')
    P = opt.P;
else
    P = size(q,1);
end
if isfield(opt,'N')
    N = opt.N;
else
    N = size(rvec_in,1)/P;
end
if isfield(opt,'M')
    M = opt.M;
else
    M = size(rvec_out,1)/P;
end
delta = opt.delta;
iterate_RFD = opt.iterate_RFD;

if P ~= size(q,1)
    error('compute_RFD: opt.P inconsistent with q');
end

[rinDown,routDown,~] = updateNodes(rvec_in,rvec_out,P,N,M,-delta/2*RandVel);
[rinUp,routUp,transVec] = updateNodes(rvec_in,rvec_out,P,N,M,delta/2*RandVel);
qUp = q+transVec;
qDown = q-transVec;

if iterate_RFD
    if isempty(precond) || ~all(isfield(precond,{'Y','UU','LL','Kin','Tblock'}))
        error('compute_RFD: precond.{Y,UU,LL,Kin,Tblock} required for iterate_RFD');
    end
    [U2, ~, ~] = solve_brownian_mobility(qUp,rinUp,routUp, ...
        precond.Y,precond.UU,precond.LL,precond.Kin,precond.Tblock, ...
        RandVel,0,opt);
    [U1, ~, ~] = solve_brownian_mobility(qDown,rinDown,routDown, ...
        precond.Y,precond.UU,precond.LL,precond.Kin,precond.Tblock, ...
        RandVel,0,opt);
else
    if isempty(precond) || ~isfield(precond,'Ktot') || ~isfield(precond,'Btot')
        % Compute Ktot/Btot if not provided
        [~,~,~,Kin,B] = oneBodyPrecondMobDLP(rvec_in(1:N,:), rvec_out(1:M,:), q(1,:));
        precond.Ktot = kron(eye(P),Kin);
        precond.Btot = kron(eye(P),B);
    end
    S_up = stokes_SLP_mat(rinUp,routUp);
    [Ytot,Utot] = getPseudoFactors(S_up,1e-13,0);
    Rup = precond.Ktot'*Ytot*(Utot'*precond.Btot);
    Mup = Rup\eye(6*P);

    S_down = stokes_SLP_mat(rinDown,routDown);
    [Ytot,Utot] = getPseudoFactors(S_down,1e-13,0);
    Rdown = precond.Ktot'*Ytot*(Utot'*precond.Btot);
    Mdown = Rdown\eye(6*P);

    U2 = Mup*RandVel;
    U1 = Mdown*RandVel;
end

RFD = (U2-U1)/delta;
end
