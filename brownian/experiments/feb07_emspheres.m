clear;
close all

solve_dense = 0;
iterate_RFD = 1; % or solve for drift term densely

rng(5);

% We test Brownian MFS by time-stepping systems at equilibrium with P = 1 or P = 2 particles 
% using Euler-Maryama and random finite differences

%% Param selection 
tsteps = 1e5; % number of time steps
tsteps = 1e4;
%tsteps = 1e3; 
%tsteps = 1e2; 
% tsteps = 1000;
dt = 1e-1; %time step size
%dt = 1e-2;

%delta = 1e-3; % parameter in RFD scheme
%delta = 5e-1;
delta = 1e-3; %pretty good with dt = 1e-2;
%delta = 5e-2;
%delta = 1; 

%Spring model for two particles
k = 10; 
l = 4; % spring relaxation length
%l = 5; 
%l = 6;
Ufunc = @(x1,x2) k*0.5*(norm(x1-x2)-l).^2;
Ffunc = @(x1,x2) k*(norm(x1-x2)-l)*(x1-x2)./norm(x1-x2);

gmres_tol = 1e-5;
tol_lanczos = 1e-5;

maxit = 100; % 

%% Set initial particle coordinates
P=2;
if P==1
    q = [0 0 0];
else
    delta = 3; %initial particle particle distance
    warning('Choose according to equilibrium distribution')
    q = [0 0 0; 2+delta 0 0];
end

%% Set resolution and discretize
%Rp = 0.33;
Rp = 0.15; %Proxy sphere radius -- very coarse resolution
N = 50;  % Number of proxy sources

% N = 100; 
% Rp = 0.30; 

% 
%Rp = 0.10;´
%N = 20; 

% Rp = 0.63; %fine resolution
% N = 700;

a = 2; %Determines oversampling factor for the collocation points
[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a); %Assign source and collocation points

opt.dt = dt;
opt.gmres_tol = gmres_tol; 
opt.gmres_verbose = 0; 
opt.delta = delta;
opt.iterate_RFD = iterate_RFD;

N = size(rvec_in,1)/P;
M = size(rvec_out,1)/P;

%% One-body preconditioning will be the same for everyone -- precompute!
% if solve_dense
%     [Y,UU,LL,Kin,B] = oneBodyPrecondMob(rvec_in(1:N,:),rvec_out(1:M,:),q(1,:)); 
%     [Ys,UUs] = oneBodyPrecondRes(rvec_in(1:N,:),rvec_out(1:M,:)); 
% end 
    [Y,UU,LL,Kin,B] = oneBodyPrecondMobDLP(rvec_in(1:N,:),rvec_out(1:M,:),q(1,:));  
   % [Ys,UUs] = oneBodyPrecondResDLP(rvec_in(1:N,:),rvec_out(1:M,:),rvec_out(1:M,:)-q(1,:)); 
%end

%for debugging to check the action of the square root of the mobility
%matrix.
%if P == 2 && solve_dense
    Ktot = blkdiag(Kin,Kin);
    Btot = blkdiag(B,B); 
%end

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
T = stokes_DLP_mat(rvec_out,rvec_in,n);




%% Prepare for preconditioned Lanczos using block diagonal preconditioning
Ti = stokes_DLP_mat(rvec_out(1:M,:),rvec_in(1:N,:),rvec_out(1:M,:)); 
Si = stokes_SLP_mat(rvec_in(1:N,:), rvec_out(1:M,:));
Ai = Ti*Si;
Ci = (Ai+Ai')/2;
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
PP = kron(eye(P),Pi);
PPplus = kron(eye(P),Pi_plus);
            

%% Loop in time
d = zeros(tsteps,1);
iter_hist = zeros(tsteps,1); 
qhist = cell(tsteps,1);

for i = 1:tsteps
    i
    
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
        precond.Ktot = Ktot;
        precond.Btot = Btot;
        opt.P = P;
        opt.N = N;
        opt.M = M;
        [RFD,~,~] = compute_RFD(q,rvec_in,rvec_out,RandVel,opt,precond);
    else
        RFD = zeros(6*P,1);
    end

    %Solve with fluctuating velocity field in rhs
    % W = randn(3*M*P,1);
    % if P==1
    %     [U, iters, lambda_norm] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,Fvec, 1, opt,sqrtA,W);
    % 
    %     % Debugging: these should be the same!
    %     %U2 = sqrt(2/dt)*(R\(Kin'*Ys*UUs*sqrtA*W));
    %     %U = sqrt(2/dt)*L*randn(6,1);
    %     %U = sqrt(2 * dt / (6*pi)) * randn(6,1);
    % else
        %% Solve with dense mobility matrix M and sqrt. The naive and dense way of solving the problem.
        G = stokes_SLP_mat(rvec_in,rvec_out);

        if solve_dense
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

            %% Solve instead mobility problem iteratively

            % Prepare for Lanczos
            A = T*G;
            C = M/(4*pi)*(A+A')/2; % correlation matix for fluctuating velocity field
            dW = randn(3*N*P,1); % white noise

            % Some debug
            %sqrtA = chol(Asym); %results in very large errors!
            
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
           % 
           %  %same thing
           %  sqrtA = Vs(:,ind)*diag(sqrt(ds(ind)))*Vs(:,ind)';
            
            % Check that we compute the right RBM!
            % [Ytot,Utot] = getPseudoFactors(A,1e-13,0);
            % Ainv = Ytot*Utot';
            % R = M/(4*pi)*Ktot'*Ainv*Ktot;
            % R2 = -Ktot'*Ainv*Tdiag*Btot;
            % [Ytot,Utot] = getPseudoFactors(S,1e-13,0);
            % R3 = Ktot'*Ytot*(Utot'*Btot);

            % Udense = (R3\eye(P*6))*(Fvec+Ktot'*Ainv*sqrtA*sqrt(2/opt.dt)*dW);
            %noise = sqrtA*dW;

            %% Set up Lanczos with preconditioning (necessary for convergence!)
            B_precond = PP*C*PP';     
            [noise,iter_errv5,iters] = KrylovSqrtMsing(B_precond,dW,tol_lanczos,maxit,PPplus);
            iter_hist(i) = iters;
    %% Solve mobility problem
           % [U, iters, lambda_norm] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,Tblock,Fvec,opt, noise);
            [U, iters, lambda_norm] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,Tblock,Fvec,opt, noise);
        end
    %end

    %Update particle coordinates (take time step) (assuming scaling such
    %that kT = 1).
    
    Unew = dt*(U+RFD);
    [rvec_in,rvec_out,transVec] = updateNodes(rvec_in,rvec_out,P,N,M,Unew); %Rotation to be taken into account? 
    % if norm(Unew,inf)>1
    %     disp('large change')
    % end
    q = q+transVec;
    if PP == 1
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



% function [rinUp,routUp,transVec] = updateNodes(rin,rout,P,N,M,vel)
% 
% transVec = zeros(P,3);
% rinUp = zeros(N*P,3);
% routUp = zeros(M*P,3);
% 
% for k = 1:P
%     velTrans = vel(6*(k-1)+1:6*(k-1)+3)';
%     transVec(k,:) = velTrans;
%     rinUp(N*(k-1)+1:k*N,:) = rin(N*(k-1)+1:k*N,:)+velTrans;
%     routUp(M*(k-1)+1:k*M,:) = rout(M*(k-1)+1:k*M,:)+velTrans;
% end
% 
% 
% end

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
