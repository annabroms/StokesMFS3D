clear; close all; 
%Builds correlation matrix for <UU'> when the noise is given by the
% action of chol(S_RPY). Compares either against the known mobility matrix,
% when this is known for a single body, or against the resulting mobility
% matrix obtained with a dense MFS solve. Conclusion: the error has a
% strong dependence on the regularisation parameter in the RPY tensor and
% two error contributions must be balanced: regularisation and discretisation. 
% Regardless, the error for the best case scenario is large!

% set particle coordinates
P=1;
delta = 5; 
if P==1
    q = [0 0 0];
else
    delta = 3; %initial particle particle distance
    q = [0 0 0; 2+delta 0 0];
end

% Set resolution
Rp = 0.63; %fine resolution
N = 700;

% Rp = 0.15; %Proxy sphere radius -- very coarse resolution
% N = 50;  % Number of proxy sources
% % 
Rp = 0.30;
N = 100; 

a = 1.2; %Determines oversampling factor for the collocation points
[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a); %Assign source and collocation points
N = size(rvec_in,1)/P;
MM = size(rvec_out,1)/P; 

B = getKmat(rvec_out,q);
K = getKmat(rvec_in,q);

G = generate_stokes_mat(rvec_in, rvec_out);
tol = 1e-13; 
visualise = 0; 
[Y,U]  = getPseudoFactors(G,tol,visualise);

R = K'*Y*(U'*B);
M = R\eye(P*6); 
M_true = diag([1/6/pi*ones(1,3) 1/8/pi*ones(1,3)]);
M-M_true

avec = logspace(-2,0);
for i = 1:length(avec)
    i
    S_RPY = generate_RPY_matrix(rvec_out,avec(i));
    S = generate_stokes_mat(rvec_out, rvec_out);

    % figure()
    % imagesc(log10(abs(S-S_RPY)))
    % colorbar
    
    KGplus = K'*Y*U';

    %Build correlation for the rigid body motion vectors.
    M_corr = M*KGplus*S_RPY*KGplus'*M';

    err(i) = norm(M-M_corr)/norm(M); 
end

%%
figure()
loglog(avec,err,'.-')
ylabel('$\|M-UU^T\|_2/\|M\|_2$','interpreter','latex')
xlabel('Regularisaiton parameter a','interpreter','latex')
