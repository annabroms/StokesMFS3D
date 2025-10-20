function [Rinv,Zi,Yi,db] = get_long_range_precond(q,rin,rout,RR,opt,V,U)
%GET_LONG_RANGE_PRECOND  Construct coarse-space projection matrices for long-range preconditioning.
%
%   [Rinv, Zi,Yi,db] = GET_LONG_RANGE_PRECOND(q, rin, rout, opt)
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
%     Zi       - (3M×k) matrix mapping coarse coefficients to proxy source strengths for a single body (diagonal block in Z).
%     Yi       - (3N×k) matrix whose transpose, Y', maps to coarse velocity space for a single body (diagonal block in Y).
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

if lr == 1
    %Use matrices Zi and Yi constructed with ones matrices
    Zi = repmat(eye(3),N,1);
    Yi = repmat(eye(3),M,1);
    db = 3; %dimension per body of the coarse source    

elseif lr == 2

    if opt.ellipsoid
        Zi = getKmat((RR{1}'*rin(1:N,:)')',q(1,:));
        Yi = getKmat((RR{1}'*rout(1:M,:)')',q(1,:));
    else
        Zi = getKmat(rin(1:N,:),q(1,:));
        Yi = getKmat(rout(1:M,:),q(1,:));
    end

    db = 6; %dimension per body of the coarse source
else
    lmax = lr-2;
    Zi = V(:,1:lmax);
    Yi = U(:,1:lmax);
    db = lmax; 
end

u = zeros(3*M*P,P*db);
for k = 1:P
        if opt.ellipsoid && (lr>1)
            Rk = RR{k};

            if lr == 2
                %block = getKmat(rin(N*(k-1)+1:k*N,:),q(k,:));
                %same thing
                block = kron(eye(N),Rk)*Zi*[Rk' zeros(3); zeros(3) Rk'];
            else
                block = reshape(rotate_vector(Zi,Rk),3*N,db);
                % this can be compared to the block Zi obtained from this
                % specific particle

                %debug mode: 
                %Gk = generate_stokes_mat(rin((k-1)*N+1:k*N,:),rout((k-1)*M+1:k*M,:));
                %do svd
                %[Uk,Sk,Vk] = svd(Gk); 
                %block_k = Vk(:,1:db); %same thing as block

                
            end

        else
            block = Zi;
        end
   
    for i = 1:db
        u(:,(k-1)*db+i) = getFlow(block(:,i),rin(N*(k-1)+1:k*N,:),rout,opt);
    end
end

for k = 1:P
    
    if lr == 2
        if opt.ellipsoid 
            Rk = RR{k};
           
            %Yik = getKmat(rout(M*(k-1)+1:k*M,:),q(k,:)); % same thing
            Yik = kron(eye(M),Rk)*Yi*[Rk' zeros(3); zeros(3) Rk'];
            R((k-1)*db+1:k*db,:) = Yik'*u((k-1)*3*M+1:k*3*M,:);
        else
            R((k-1)*db+1:k*db,:) = Yi'*u((k-1)*3*M+1:k*3*M,:);
        end

    elseif lr>2
        if opt.ellipsoid
            Rk = RR{k};
            Yik = reshape(rotate_vector(Yi,Rk),3*M,db);
            R((k-1)*db+1:k*db,:) = Yik'*u((k-1)*3*M+1:k*3*M,:);
        else
            R((k-1)*db+1:k*db,:) = Yi'*u((k-1)*3*M+1:k*3*M,:);
        end

    else
        R((k-1)*db+1:k*db,:) = Yi'*u((k-1)*3*M+1:k*3*M,:);
    end

    
end

% Faster computation using FMM??

Rinv = (R\eye(db*P));
   
    
end

