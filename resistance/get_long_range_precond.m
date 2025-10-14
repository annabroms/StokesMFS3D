function [Rinv,Z,Y] = get_long_range_precond(q,rin,rout,opt)
%GET_LONG_RANGE_PRECOND  Construct coarse-space projection matrices for long-range preconditioning.
%
%   [Rinv, Z, Y] = GET_LONG_RANGE_PRECOND(q, rin, rout, opt)
%
%   Constructs the matrices used in the long-range preconditioner:
%   the coarse-to-fine mappings AN and AM, and the inverse coarse interaction matrix Rinv.
%   These are used to define the projectors:
%
%       P = I - G * Z * Rinv * Y'     (left-side projection, see applyPmat)
%       Q = I - Z * Rinv * Y' * G     (right-side projection, see applyQmat)
%
%   The choice of coarse basis (e.g. translation-only vs. rigid-body modes) 
%   is set via the option opt.lr.
%
%   INPUTS:
%     q        - P×3 array of particle centers.
%     rin      - PN×3 array of proxy/source points (used here to determine for coarse basis flow).
%     rout     - PM×3 array of surface/target points (used to evaluate coarse flow).
%     opt      - Struct with fields:
%                  • lr: 1 = translation-only modes (3 per body),
%                        2 = full rigid-body modes (6 per body),
%                  • N: number of proxy points per body,
%                  • M: number of collocation points per body,
%                  • other kernel options (such as fmm flat passed to getFlow).
%
%   OUTPUTS:
%     Rinv     - (3Pk×3Pk) inverse of coarse interaction matrix R = AM'*U,
%                where U is the coarse-flow matrix at the PM targets from
%                the Pk coarse basis functions
%     Z       - (3PM×3Pk) matrix mapping coarse coefficients to proxy source strengths.
%     Y       - (3PN×3Pk) matrix whose transpose, Y', maps to coarse velocity space.
%
%   NOTES:
%     - The function assumes block structure per body and constructs
%       block-diagonal matrices Z and Y accordingly.
%     - Coarse-flow matrix U is built using calls to getFlow with basis sources.
%     - Currently uses dense matrix assembly; consider FMM for acceleration.
%
%   See also: applyPmat, applyQmat, getFlow, getKmat
%
% Anna Broms, Oct 14, 2025

lr = opt.lr;
M = opt.M; 
N = opt.N; 
P = size(q,1);

%Use different matrices Z and Y using ones matrices instead!
if lr == 1
    for k = 1:P % not the fastest way to generate this...
        %ones in the x,y,z direction, both for Z and Y
        blocksN{k} = repmat(eye(3),N,1);
        blocksM{k} = repmat(eye(3),M,1);
       
    end  
    db = 3; %dimension per body of the coarse source    

elseif lr == 2
    for k = 1:P
        blocksN{k} = getKmat(rin(N*(k-1)+1:k*N,:),q(k,:)); % apply rotations intead?
        blocksM{k} = getKmat(rout(M*(k-1)+1:k*M,:),q(k,:));
    end
    db = 6; %dimension per body of the coarse source
else
    error("not yet implemented")
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if lr ==1
    block = repmat(eye(3),N,1);
    s=3;
else
    s=6;
end

for k = 1:P
    if lr == 1
        blocks{k} = block;
    elseif lr == 2
        block = getKmat(rin(N*(k-1)+1:k*N,:),q(k,:));
        blocks{k} = block;
    end
    for i = 1:s
        u(:,(k-1)*s+i) = getFlow(block(:,i),rin(N*(k-1)+1:k*N,:),rout,opt);
    end
end

Z = blkdiag(blocksN{:});
Y = blkdiag(blocksM{:});

R = Y'*u;

% Faster computation using FMM??

Rinv = (R\eye(db*P));
   
    
end

