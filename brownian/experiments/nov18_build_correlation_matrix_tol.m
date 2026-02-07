clear; close all;
% Build correlation matrix for <UU'> when noise is given by fluctuating
% sources with correlation TG. This script checks the error in the
% correlation matrix <UU^T> relative to the known mobility matrix.
%
% Overview:
%   1) Build dense operators (G, T) and mobility matrix M.
%   2) Form TG = T*G and its symmetric part C = (TG+TG')/2.
%   3) For a sweep of truncation tolerances, build pseudoinverses and
%      compare the resulting correlation matrices against M_true.
%   4) Report symmetry error in the pseudoinverse and plot convergence as
%      function of truncation tolerance
%
% Notes:
%   Script is intended for verification, based on dense operators.

%% Set particle coordinates and resolution
P=2;
delta = 1; 
if P==1
    q = [0 0 0];
else
    delta = 3; %initial particle particle distance
    q = [0 0 0; 2+delta 0 0];
end

% Set resolution
% Rp = 0.63; %fine resolution
% N = 700;

Rp = 0.15; %Proxy sphere radius -- very coarse resolution
N = 50;  % Number of proxy sources
% % 
Rp = 0.30;
N = 100; 

a = 1.2; %Determines oversampling factor for the collocation points

%% Discretise particles
[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a); %Assign source and collocation points
N = size(rvec_in,1)/P;
MM = size(rvec_out,1)/P; 

%% Get all matrices needed to solve the system densely
B = getKmat(rvec_out,q);
K = getKmat(rvec_in,q);
n = repmat(rvec_out(1:MM,:),P,1);

% if P == 1
%     %Tt = getTractionMat(rvec_in,rvec_out,n); %transpose of T
%     T = stokes_DLP_mat(rvec_out,rvec_in,n);
% else
%     %Ttblock = getTractionMat(rvec_in(1:N,:),rvec_out(1:MM,:),rvec_out(1:MM,:));
%     %T = kron(eye(P),Tblock);
%     T2 = getTraction(rvec_in,rvec_out,n);
% end

T = stokes_DLP_mat(rvec_out,rvec_in,n);
G = stokes_SLP_mat(rvec_in, rvec_out);
visualise = 1; 
tol = 1e-13;

A = 4*pi/MM*T*G;

%% Determine factorisations and flavors of correlation matrices
[YG,UG]  = getPseudoFactors(G,tol,visualise);

% Determine mobility matrix as in standard MFS
R = K'*YG*(UG'*B);
M = R\eye(P*6); 

% Reference mobility matrix
if P == 2
    M_true = M;
else
    M_true = diag([1/6/pi*ones(1,3) 1/8/pi*ones(1,3)]);
    M-M_true
end

%G (SLP) pseudoinverse
Gplus = YG*UG';

%Compare eigvals of TG to those of the RPY tensor
[VA,DA] = eig(A);
d = diag(DA);
figure()
semilogy(abs(d));
title('Eigvals of TG and RPY')

S_RPY  = generate_RPY_matrix(rvec_out,0.1);
[V_RPY,D_RPY] = eig(S_RPY); 
hold on
semilogy(diag(D_RPY));

% Look at symmetric part of TG
C = (A+A')/2;

Cinv = C\eye(size(C)); %Not well-behaced due to ill-conditioning!
% %B = chol(C); 
% R_b = (K'*Cinv*K);
% M_b = R_b\eye(6*P); 
% M_corr_c = M_b*K'*Cinv*C*Cinv'*K*M_b'; % The algebra for UU^T using this mobility matrix

% Ainv = A\eye(size(A)); % Not well-behaved due to ill-conditioning
% R_b = (K'*Ainv*K);
% M_b = R_b\eye(6*P);
% M_corr_d = M_b*K'*Ainv*C*Ainv'*K*M_b'; % UU^T with this choice: symmetric part of TG used for the noise.

% Preparation for truncating eigvals
[VS,DS] =  eig(C);
ds = diag(DS); 


tolvec = logspace(-15,-1,20);
for i = 1:length(tolvec)
    i
    tol = tolvec(i);
    [YT,UT]  = getPseudoFactors(A,tol,0);
    TGplus = YT*UT';
    AA = TGplus-TGplus';
    aa(i) = norm(AA,inf); %this is the symmetry error in the pseudoinverse
    %eTG_plus = eig(TGplus); 


    %Approach 1a based in eigendecomp
   %  ind = find(real(d)>tol);
   %  diff_set = setdiff(1:length(d),ind);
   %  d2 = 1./d; 
   %  d2(diff_set) = 0;
   % % A2 = V*diag(d2)*inv(V);
   %  TG_inv = V*diag(d2)*(V\eye(size(V)));
   %  R_b = K'*TG_inv*K;
   %  M_b = R_b\eye(6*P);
   % 
   %  M_corr1a = M_b*K'*TG_inv*A*TG_inv'*K*M_b';

    % Use eigendecomp of symmetriced matrix, truncate
    ind = find(real(ds)>tol);
    diff_set = setdiff(1:length(ds),ind);
    invd2 = 1./ds; 
    invd2(diff_set) = 0;

    d2 = ds;
    d2(diff_set) = 0;
   % A2 = V*diag(d2)*inv(V);
    TG_inv = VS*diag(invd2)*VS';
    TG_reg = VS*diag(d2)*VS';
    R_b = (K'*TGplus*K);
    M_b = R_b\eye(6*P);
    %verison 3 of the mobility matrix
    
    %Use truncated eigvals of symmetric part for the noise
    M_corr1a = M_b*K'*TGplus*TG_reg*TGplus'*K*M_b';
    erra(i) = norm(M_true-M_corr1a,2)/norm(M_true,2);

    %... or don't truncate the noise matrix
    M_corr1b = M_b*K'*TGplus*C*TGplus'*K*M_b';
    errb(i) = norm(M_true-M_corr1b,2)/norm(M_true,2);

    %Approach 1b based on Eig decomp 

    %M_corr1b = M_b*K'*TGplus*A*TGplus'*K*M_b';
   % errb(i) = norm(M_true-M_corr1b,2)/norm(M_true,2);

    R_b = K'*TG_inv*K;
    M_b = R_b\eye(6*P);
    M_corr1c = M_b*K'*TG_inv*TG_reg*TG_inv'*K*M_b';
    errc(i) = norm(M_true-M_corr1c,2)/norm(M_true,2);

   
    %Approach 2: 
    R_2 = (K'*Gplus*B);
    M_2 = R_2\eye(6*P);
    % M_corr2 = M_2*K'*TGplus*A*TGplus'*K*M_2';
    % err2(i) = norm(M_true-M_corr2,2)/norm(M_true,2);

    %Approach 3:
    % M_corr3 = M_2*K'*Gplus*S*TGplus*S'*Gplus'*K*M_2';
    % err3(i) = norm(M_true-M_corr3,2)/norm(M_true,2);

    %Approach 4 (standard solve of problem)
    M_corr4 = M_2*K'*Gplus*G*TG_inv*G'*Gplus'*K*M_2';
    err4(i) = norm(M_true-M_corr4,2)/norm(M_true,2);

    M_corr5 = M_2*K'*Gplus*G*Cinv*G'*Gplus'*K*M_2';
    err5(i) = norm(M_true-M_corr4,2)/norm(M_true,2);

    

end
%% Plot correlation errors vs truncation tolerance
figure()
loglog(tolvec,erra,'+-','DisplayName','Pseduoinv (SVD) TG solve -- eigendecomp sym(TG) noise');
hold on
loglog(tolvec,errb,'*-','DisplayName','Pseduoinv (SVD) TG solve -- TG noise (no regularisation)');
loglog(tolvec,errc,'.-','DisplayName','Pseduoinv (Eigen) sym(TS) solve -- eigendecomp sym(TS) noise');
%loglog(tolvec,err2,'*-');
%loglog(tolvec,err3,'o-');
loglog(tolvec,err4,'o-','DisplayName','Pseudoinv (SVD) G solve -- eigendecomp TG noise');
loglog(tolvec,err5,'o-','DisplayName','Pseudoinv (SVD) G solve -- eigendecomp sym(TG) noise');
axis tight
xlabel('eig/svd truncation TOL','Interpreter','latex')
ylabel('$\|M-UU^T\|_2/\|M\|_2$','interpreter','latex')
ylim([1e-11,1])
grid on
legend
set(legend,'interpreter','latex')
% Do the same with an eigenvalue factorization

figure()
loglog(tolvec,aa)
xlabel('eig/svd truncation TOL','Interpreter','latex')
ylabel('Symmetry error in SVD pseudoinverse of TS','Interpreter','latex')
axis tight
grid on 
