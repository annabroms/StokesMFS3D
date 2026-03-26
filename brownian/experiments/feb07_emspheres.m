clear;
close all

solve_dense = 0;
iterate_RFD = 0; % or solve for drift term densely
use_dense_lanczos = 0; % use dense matrices in Lanczos (old path) vs matvecs

rng(5);

% We test Brownian MFS by time-stepping systems at equilibrium with P = 1 or P = 2 particles 
% using Euler-Maryama and random finite differences

%% Param selection 
tsteps = 1e5; % number of time steps
tsteps = 1e4;
%tsteps = 1e3; 
%tsteps = 1e2; 
dt = 1e-1; %time step size
%dt = 1e-2;

%delta_rfd = 1e-2; % parameter in RFD scheme
delta_rfd = 1e-3; %pretty good with dt = 1e-2;
%delta_rfd = 5e-2;
 

%Spring model for two particles
k = 10; 
l = 4; % spring relaxation length
%l = 5; 
%l = 6;
Ufunc = @(x1,x2) k*0.5*(norm(x1-x2)-l).^2;
Ffunc = @(x1,x2) k*(norm(x1-x2)-l)*(x1-x2)./norm(x1-x2);

gmres_tol = 1e-4;
tol_lanczos = 1e-4;

maxit = 100; % 

%% Set initial particle coordinates
P=2;
if P==1
    q = [0 0 0];
else
    delta_init = 3; %initial particle particle distance
    warning('Choose according to equilibrium distribution')
    q = [0 0 0; 2+delta_init 0 0];
end

%% Set resolution and discretize
%Rp = 0.33;
Rp = 0.15; %Proxy sphere radius -- very coarse resolution
N = 50;  % Number of proxy sources

% N = 100; 
% Rp = 0.30; 

% 
% Rp = 0.10;
% N = 20; 

% Rp = 0.63; %fine resolution
% N = 700;

a = 2; %Determines oversampling factor for the collocation points
[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a); %Assign source and collocation points

opt.dt = dt;
opt.gmres_tol = gmres_tol; 
opt.gmres_verbose = 0; 
opt.delta = delta_rfd;
opt.iterate_RFD = iterate_RFD;

N = size(rvec_in,1)/P;
M = size(rvec_out,1)/P;
w = (4*pi/M)*ones(P*M,1);

%% One-body preconditioning will be the same for everyone -- precompute!
[Y,UU,LL,Kin,B] = oneBodyPrecondMobDLP(rvec_in(1:N,:),rvec_out(1:M,:),q(1,:));  

% Build block-diagonal K/B only when needed (dense solve or dense RFD)
need_dense_blocks = solve_dense || ~iterate_RFD;
if need_dense_blocks
    Ktot = kron(eye(P),Kin);
    Btot = kron(eye(P),B);
end

%% If dense solve, stuff can be precomputed
% if P == 1
%     % Get sqrt of mobility matrix for single body
%     R = Kin'*Ys*(UUs*B);
%     MM = R\eye(6);
%     L1 = chol(MM);
%     L = diag(sqrt([1/6/pi*ones(3,1);1/8/pi*ones(3,1)])); %what it should be
% 
%     % Precompute sqrt
%     %[u,U,ss_sqrt,V] = getRandomField(q,rvec_in,rvec_out,P); %old -- does
%     %not give the right correlations
%     %sqrtA = U*(ss_sqrt*V');
% end

Tblock = stokes_DLP_mat(rvec_out(1:M,:),rvec_in(1:N,:),rvec_out(1:M,:)-q(1,:));
n = repmat(rvec_out(1:M,:)-q(1,:),P,1);
if use_dense_lanczos
    T = stokes_DLP_mat(rvec_out,rvec_in,n);
end

%% Prepare for preconditioned Lanczos using block diagonal preconditioning
Ti = stokes_DLP_mat(rvec_out(1:M,:),rvec_in(1:N,:),rvec_out(1:M,:)); 
Si = stokes_SLP_mat(rvec_in(1:N,:), rvec_out(1:M,:));
Ai = Ti*Si;
Ci = M/(4*pi)*(Ai+Ai')/2;
[Ve,De] = eig(Ci);
d = diag(De);

ind = find(real(d)>1e-14);
diff_set = setdiff(1:length(d),ind);
%d2 = sqrt(d(ind));

dsqrt = 1./sqrt(d);
dsqrt(diff_set) = 0;
Pi = diag(dsqrt)*Ve';
dsqrt_plus = sqrt(d); 
dsqrt_plus(diff_set) = 0;
Pi_plus = Ve*diag(dsqrt_plus);

PPplus = @(x) apply_block_diag(x, Pi_plus, P);            

%% Loop in time
d = zeros(tsteps,1);
iter_hist = zeros(tsteps,1); 
qhist = cell(tsteps,1);

for i = 1:tsteps
    i
    
    %% Set external force
    if P == 2
        d(i) = norm(q(1,:)-q(2,:)); %compute center center distance
        Fvec = [-Ffunc(q(1,:),q(2,:))'; zeros(3,1); Ffunc(q(1,:),q(2,:))'; zeros(3,1)];
        qhist{i} = q;

        %Fvec = zeros(P*6,1); %test
    else
        Fvec = zeros(P*6,1);
        xhist(i,:) = q;
        
    end

    %% Compute random finite difference contribution
    if P == 2
        RandVel = randn(6*P,1);
        precond.Y = Y;
        precond.UU = UU;
        precond.LL = LL;
        precond.Kin = Kin;
        precond.Tblock = Tblock;
        if need_dense_blocks
            precond.Ktot = Ktot;
            precond.Btot = Btot;
        end
        opt.P = P;
        opt.N = N;
        opt.M = M;
        [RFD,~,~] = compute_RFD(q,rvec_in,rvec_out,RandVel,opt,precond);
    else
        RFD = zeros(6*P,1);
    end

   
    %% Solve with dense mobility matrix M and sqrt. The naive and dense way of solving the problem.
    if solve_dense
        G = stokes_SLP_mat(rvec_in,rvec_out);
        %build sqrt of M densely:
        [Ytot,Utot] = getPseudoFactors(G,1e-13,0);
        Rtot = Ktot'*Ytot*(Utot'*Btot);
        Mtot = Rtot\eye(6*P);
        U = Mtot*Fvec;
        try 
            R = chol(Mtot);
        catch
            disp('Not pos def')
        end
        W = sqrt(2/dt)*randn(6*P,1);
        U = U+R'*W;
    else
        %% Solve instead with fluctuating velocity field in rhs, with the
        % mobility solve iteratively

        % Prepare for Lanczos
        dW = randn(3*N*P,1); % white noise

       % Some debug
       % sqrtA = chol(Asym); %results in very large errors!
        
       %  [Vs,Ds] = eig(C);
       % % Ainv = A\eye(size(A)); 
       %  ds = diag(Ds); 
       %  tol = 1e-8; %works
       %  tol = 1e-15; %1e-12 works
       %  ind = find(real(ds)>tol);
       %  % diff_set = setdiff(1:length(ds),ind);
       %  % sqds = sqrt(ds); 
       %  % sqds(diff_set) = 0;
       %  % sqrtA = Vs*diag(sqds)*Vs';
   
        %noise = sqrtA*dW;

        %% Set up Lanczos with preconditioning (necessary for convergence!)
        if use_dense_lanczos
            G = stokes_SLP_mat(rvec_in,rvec_out);
            A = T*G;
            C = M/(4*pi)*(A+A')/2; % correlation matrix for fluctuating velocity field
            B_precond = PP*C*PP';     
           % [noise,iter_errv,iters] = lanczosSqrt(B_precond,dW,tol_lanczos,maxit);
            [noise,iter_errv,iters] = lanczosSqrt(B_precond,dW,tol_lanczos,maxit,PPplus);
           % [noise2,iter_errv,iters2] = lanczosSqrt_fast(B_precond,dW,tol_lanczos,maxit,PPplus);

            Mfun = @(x) getPrecondTG(x,P,rvec_in,rvec_out,n,w,Pi,opt);
            %[noise3,iter_errv,iters3] = lanczosSqrt_fast(Mfun,dW,tol_lanczos,maxit,PPplus);
           % disp('debug')
        else
            % Build matrix via matvec
            Mfun = @(x) getPrecondTG(x,P,rvec_in,rvec_out,n,w,Pi,opt);
            [noise,iter_errv,iters] = lanczosSqrt_fast(Mfun,dW,tol_lanczos,maxit,PPplus);
      
        end
        iter_hist(i) = iters;
      %% Solve mobility problem with modified boundary data
        [U, iters, lambda_norm] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,Tblock,Fvec,opt, noise);
    end

    %Update particle coordinates (take time step) (assuming for now scaling such
    %that kT = 1).
    
    Unew = dt*(U+RFD);
    [rvec_in,rvec_out,transVec] = updateNodes(rvec_in,rvec_out,P,N,M,Unew); %Rotation to be taken into account? 

    q = q+transVec;
    if P == 1
        d(i) = norm(q)^2; %store squared displacement
    end    
    
end

figure()
scatter3(rvec_in(:,1),rvec_in(:,2),rvec_in(:,3),'magenta','filled')
hold on
scatter3(rvec_out(:,1),rvec_out(:,2),rvec_out(:,3),'blue','filled')
axis equal

%% Visualise mean square displacement
if P == 1
    %D = mean(d)/(2*dt*i);
    Dhist = cumsum(d)'./(1:i);    
    figure()
    plot(dt*(1:i),1/6/pi*ones(1,i),'k--');
    hold on
    plot(dt*(1:i-1),Dhist(2:end)./(2*dt*(1:i-1)),'b-')
    [msd,lags] = compute_msd(xhist);

    %%
    figure()
    subplot(1,2,1)
    Nshow = 10;
    plot(dt*lags,msd');
    subplot(1,2,2);
    plot(dt*lags,msd'./(6*(dt*lags)))  
    hold on
    plot(dt*lags,1/6/pi*ones(size(lags)),'k--');
    ylim([0,0.2])

    %% Plot MSD
    figure;
    loglog(dt * lags, msd, 'b', 'LineWidth', 2);
    hold on;
    loglog(dt * lags, 1/pi * dt * lags, 'k--', 'LineWidth', 2);
    xlabel('Lag time \tau'); ylabel('MSD(\tau)');
    legend('Simulated MSD', 'Expected 6D\tau');
    title('Brownian motion of a sphere: MSD test');
    grid on;
else  %P == 2 -- Draw marginal pdf
    figure()
    histogram(d(1:i),50,'Normalization','pdf')
    hold on
    k = 10;
    p = @(r) 4*pi*r.^2/((2*pi/k).^(3/2)).*exp(-k*(r-l).^2/2); %for the potential spring
    p2 = @(r) r.^2.*exp(-k*(r-l).^2/sqrt(2)); %for the potential spring
    rmax = 7;
    Q = integral(p,1,rmax);
    pn = @(r) p(r)/Q;
    Q2 = integral(p2,1,rmax);
    pn2 = @(r) p2(r)/Q2;
    r = linspace(2,rmax,200);
    plot(r,pn(r),'k--','LineWidth',2);
    xlabel('particle-particle distance $d$','Interpreter','latex')
    ylabel('$p(d)$','Interpreter','latex')
   % plot(r,pn2(r),'r--')
end



function [msd,lags] = compute_msd(X)
    % X: N x d array of positions over time
    N = size(X, 1);
    maxLag = floor(N/2);
    %maxLag = 3*N/4;
   % maxLag = N; 
    lags = 1:100:floor(maxLag);
    nlags = length(lags);
    msd = zeros(nlags, 1);
    
    for i = 1:nlags
        n = lags(i);
        diffs = X(1 + n:N, :) - X(1:N - n, :);
        sq_displacements = sum(diffs.^2, 2);
        msd(i) = mean(sq_displacements);
    end
end
    
