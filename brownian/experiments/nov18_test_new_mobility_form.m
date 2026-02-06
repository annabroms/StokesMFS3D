close all; clear;

%Can we solve a resistance problem by first multiplying through with the
%DLP operator Tt? The experiment in this scripts suggests that this is
%possible. We quantify both the error in the resistance matrix and the
%force error as some random rigid body velocity determines the boundary
%data. The matrix Tt is supposed to be a dense matrix, but here we also
%played with Tt block diagonal. 

%% Set initial particle coordinates
P=2;
delta = 1; 
if P==1
    q = [0 0 0];
else
    q = [0 0 0; 2+delta 0 0];
end
fprintf('Solve dense resistance problems for P = %d particle(s), delta = %.3g apart\n', P, delta);

%% Discretize the particle(s)
%Not surprisingly, the accuracy will depend on resolution.
%Rp = 0.33;
Rp = 0.15; %Proxy sphere radius -- very coarse resolution
N = 50;  % Number of proxy sources
% 
% % 
% Rp = 0.10;
% N = 20; 

% Rp = 0.63; %fine resolution
% N = 700;

a = 3; %Determines oversampling factor for the collocation points
[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a); %Assign source and collocation points
N = size(rvec_in,1)/P;
M = size(rvec_out,1)/P; 

%% Get all the necessary discretised operators
%rvec_in = rvec_in(1,:);
n = rvec_out - kron(q,ones(M,1));
if P == 1
    Tt = getTractionMat(rvec_in,rvec_out,n);
    T2 = Tt;
else
    %Tt block diagonal here
    Tblock = getTractionMat(rvec_in(1:N,:),rvec_out(1:M,:),rvec_out(1:M,:)-q(1,:));
    Tt = kron(eye(P),Tblock);
    T2 = getTractionMat(rvec_in,rvec_out,n);
end
G = stokes_SLP_mat(rvec_in, rvec_out); %to match naming convention of document.

visualise = 1; 
tol = 1e-10; %SVD truncation level
tol = 1e-15;

A = Tt'*G;
A2 = T2'*G;

[Y_G,U_G]  = getPseudoFactors(G,tol,visualise);
[Y_TG,U_TG]  = getPseudoFactors(A,tol,visualise);
tol = 1e-16; 
[Y_T2G,U_T2G]  = getPseudoFactors(A2,tol,visualise);
Kout = getKmat(rvec_out,q);
Kin = getKmat(rvec_in,q);

%% The near SPD matrix A = TG is rank deficient 
% and we therefore look only at a representation where the eigenvalues above 
% TOL are truncated.
[V,D] = eig(A);
d = diag(D);
semilogy(abs(d));
tol = 1e-9;
dprime = d(real(d)>tol);
ind = find(real(d)>tol);
diff_set = setdiff(1:length(d),ind);
d2 = 1./d;
%d2 = d; 
d2(diff_set) = 0;
inv_V = V\eye(size(V)); %stupid way of attempting this result
A3 = V*diag(d2)*inv_V; 

%% Compare three different ways to compute the forces and look at errors
Urand = rand(6*P,1); %random translational and angular velocities
F_orig = Kin'*Y_G*(U_G'*Kout*Urand);
F_traction1 = M*Kin'*Y_TG*(U_TG'*Kin*Urand)/4/pi;
F_traction2 = Kin'*Y_TG*(U_TG'*Tt'*Kout*Urand);
F_traction3 = M*Kin'*A2*Kin*Urand/4/pi; 

fprintf('\nForce Error (inf-norm, relative)\n');
fprintf('  SVD(G) (standard MFS) vs. SVD(TG): % .3e\n', norm(F_orig-F_traction1,inf)/norm(F_orig,inf));
fprintf('  SVD(G) (standard MFS) vs. SVD(TG)*T*Kout:     % .3e\n', norm(F_orig-F_traction2,inf)/norm(F_orig,inf));
%norm(F_orig-F_traction3,inf)/norm(F_orig,inf)

% What does the resistance matrix look like?
R_orig = Kin'*Y_G*(U_G'*Kout);
R_trac1 = M*Kin'*Y_TG*(U_TG'*Kin)/4/pi;
R_trac2 = Kin'*Y_TG*(U_TG'*Tt'*Kout);
R_trac3 = M*Kin'*A3*Kin/4/pi;
R_trac4 = Kin'*Y_T2G*(U_T2G'*Tt'*Kout);
R_trac5 = M*Kin'*Y_T2G*(U_T2G'*Kin)/4/pi;

fprintf('\nResistance Matrix Error (2-norm, relative)\n');
fprintf('  SVD(G) (standard MFS) vs. SVD(TG):        % .3e\n', norm(R_orig-R_trac1,2)/norm(R_orig,2));
fprintf('  SVD(G) (standard MFS) vs. SVD(TG)*T*Kout:            % .3e\n', norm(R_orig-R_trac2,2)/norm(R_orig,2));
fprintf('  SVD(G) (standard MFS) vs. trunc eig(TG):        % .3e\n', norm(R_orig-R_trac3,2)/norm(R_orig,2));
fprintf('  SVD(G) (standard MFS) vs. SVD(T2*G)*T*Kout: % .3e\n', norm(R_orig-R_trac4,2)/norm(R_orig,2));
fprintf('  SVD(G) (standard MFS) vs. SVD(T2*G):      % .3e\n', norm(R_orig-R_trac5,2)/norm(R_orig,2));

%% Look for accuracy in flow field: evaluate on surface in a new set of points
if P == 1
    lambda = M*Y_TG*(U_TG'*Kin*Urand)/4/pi;
    % Get new nodes for evaluating velocity residuals
    b = ellipsoid_param(1,1,1);   % baseline object at the origin, aligned
    b = setupsurfquad(b,[46,55]); %returns discretisation of sphere with finer grid
    rcheck = b.x';
    %rcheck = rvec_out; 
    S = stokes_SLP_mat(rvec_in, rcheck);
    res = S*lambda; 
    Kcheck = getKmat(rcheck,[0 0 0]); 
    res2 = Kcheck*Urand;
    
    diff_vec = res-res2;
    fprintf('\nSurface Velocity Residual (inf-norm, relative)\n');
    fprintf('  TG vs. rigid body: % .3e\n', norm(diff_vec,inf)./norm(res2,inf));
end
