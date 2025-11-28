close all; clear;

%Can we solve a resistance problem by first multiplying through with the
%DLP operator T? The experiment in this scripts suggests that this is
%possible. We quantify both the error in the resistance matrix and the
%force error as some random rigid body velocity determines the boundary
%data.

%% Set initial particle coordinates
P=2;
delta = 5; 
if P==1
    q = [0 0 0];
else
    delta = 3; %initial particle particle distance
    q = [0 0 0; 2+delta 0 0];
end

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
n = repmat(rvec_out(1:M,:),P,1);
if P == 1
    T = getTraction(rvec_in,rvec_out,n);
else
    %Remember T should be block diagonal
    Tblock = getTraction(rvec_in(1:N,:),rvec_out(1:M,:),rvec_out(1:M,:));
    T = kron(eye(P),Tblock);
end
S = generate_stokes_mat(rvec_in, rvec_out);

visualise = 1; 
tol = 1e-10;
tol = 1e-15;

A = -T*S;

[Y,U]  = getPseudoFactors(S,tol,visualise);
[YT,UT]  = getPseudoFactors(A,tol,visualise);
Kout = getKmat(rvec_out,q);
Kin = getKmat(rvec_in,q);

%% The near SPD matrix A = DS is rank deficient 
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
A2 = V*diag(d2)*inv_V; 

%% Compare three different ways to compute the forces and look at errors
Urand = rand(6*P,1); %random translational and angular velocities
F_orig = Kin'*Y*(U'*Kout*Urand);
F_traction1 = M*Kin'*YT*(UT'*Kin*Urand)/4/pi;
F_traction2 = -Kin'*YT*(UT'*T*Kout*Urand);
F_traction3 = M*Kin'*A2*Kin*Urand/4/pi; 

norm(F_orig-F_traction1,inf)/norm(F_orig,inf)
norm(F_orig-F_traction2,inf)/norm(F_orig,inf)
%norm(F_orig-F_traction3,inf)/norm(F_orig,inf)

% What does the resistance matrix look like?
R_orig = Kin'*Y*(U'*Kout);
R_trac1 = M*Kin'*YT*(UT'*Kin)/4/pi;
R_trac2 = -Kin'*YT*(UT'*T*Kout);
%R_trac3 = M*Kin'*A2*Kin/4/pi;

norm(R_orig-R_trac1,2)/norm(R_orig,2)
norm(R_orig-R_trac2,2)/norm(R_orig,2)
%norm(R_orig-R_trac3,2)/norm(R_orig,2)

%% Look for accuracy in flow field: evaluate on surface in a new set of points
if P == 1
    lambda = M*YT*(UT'*Kin*Urand)/4/pi;
    % Get new nodes for evaluating velocity residuals
    b = ellipsoid_param(1,1,1);   % baseline object at the origin, aligned
    b = setupsurfquad(b,[46,55]); %returns discretisation of sphere with finer grid
    rcheck = b.x';
    %rcheck = rvec_out; 
    S = generate_stokes_mat(rvec_in, rcheck);
    res = S*lambda; 
    Kcheck = getKmat(rcheck,[0 0 0]); 
    res2 = Kcheck*Urand;
    
    diff_vec = res-res2;
    norm(diff_vec,inf)./norm(res2,inf)
end