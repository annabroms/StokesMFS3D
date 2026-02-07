clear;
close all

solve_dense = 0;
iterate_RFD = 1; % or solve densely

rng(5);

% We test Brownian MFS by time-stepping systems at equilibrium with P = 1 or P = 2 particles 
% using Euler-Maryama and random finite differences

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

Tblock = getTraction(rvec_in(1:N,:),rvec_out(1:M,:),rvec_out(1:M,:)-q(1,:));
Tdiag = kron(eye(P),Tblock);
n = repmat(rvec_out(1:M,:)-q(1,:),P,1);
T = getTraction(rvec_in,rvec_out,n);


%% Param selection 
%loop in time 
tsteps = 1e5;
tsteps = 1e4;
%tsteps = 1e3; 
%tsteps = 1e2; 
% tsteps = 1000;
dt = 1e-1;
%dt = 1e-2;
 
%dt = 1e-6;
%delta = 1e-3; % parameter in RFD scheme
%delta = 5e-1;
delta = 1e-3; %pretty good with dt = 1e-2;
%delta = 5e-2;
%delta = 1; 
qhist = cell(tsteps,1);

%Spring model for two particles
k = 10; 
l = 4; % spring relaxation length
%l = 5; 
%l = 6;
Ufunc = @(x1,x2) k*0.5*(norm(x1-x2)-l).^2;
Ffunc = @(x1,x2) k*(norm(x1-x2)-l)*(x1-x2)./norm(x1-x2);

opt.dt = dt;
opt.gmres_tol = 1e-5;

tol_krylov = 1e-5;
maxit = 100; 

%% Prepare for preconditioned Lanczos using block diagonal preconditioning
Tti = getTractionMat(rvec_in(1:N,:),rvec_out(1:M,:),rvec_out(1:M,:)); 
Si = stokes_SLP_mat(rvec_in(1:N,:), rvec_out(1:M,:));
Ai = Tti'*Si;
Bi = (Ai+Ai')/2;
[Ve,De] = eig(Bi);
d = diag(De);

ind = find(real(d)>1e-14);
diff_set = setdiff(1:length(d),ind);
%d2 = sqrt(d(ind));

dsqrt = 1./sqrt(d);
dsqrt(diff_set) = 0;
Ci = diag(dsqrt)*Ve';

dsqrt_plus = sqrt(d); 
dsqrt_plus(diff_set) = 0;
Ci_plus = Ve*diag(dsqrt_plus);
C = kron(eye(P),Ci);
Cplus = kron(eye(P),Ci_plus);
            

%% Loop in time
d = zeros(tsteps,1);
iter_hist = zeros(tsteps,1); 


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

    %Compute random finite difference contribution
    if P == 2
        RandVel = randn(6*P,1);
      %  [U1, ~, ~] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,RandVel,0, opt);
             
        [rinDown,routDown,~] = updateNodes(rvec_in,rvec_out,P,N,M,-delta/2*RandVel);
        [rinUp,routUp,transVec] = updateNodes(rvec_in,rvec_out,P,N,M,delta/2*RandVel);

        if iterate_RFD
            qUp = q+transVec;
            qDown = q-transVec;
            [U2, ~, ~] = solve_brownian_mobility(qUp,rinUp,routUp,Y,UU,LL,Kin,Tblock,RandVel,0, opt);
            [U1, ~, ~] = solve_brownian_mobility(qDown,rinDown,routDown,Y,UU,LL,Kin,Tblock,RandVel,0, opt);
       
            % Solve with the standard mobility formulation instead, for debug
            % purposes. Gives the same result.
          % [U2, ~, ~] = solve_mobility(qUp,rinUp,routUp,RandVel,opt);
          % [U1, ~, ~] = solve_mobility(qDown,rinDown,routDown,RandVel,opt)
        
        else

            %Build densely
            S_up = stokes_SLP_mat(rinUp,routUp);
            [Ytot,Utot] = getPseudoFactors(S_up,1e-13,0);
            Rup = Ktot'*Ytot*(Utot'*Btot);
            Mup = Rup\eye(6*P);

            S_down = stokes_SLP_mat(rinDown,routDown);
            [Ytot,Utot] = getPseudoFactors(S_down,1e-13,0);
            Rdown = Ktot'*Ytot*(Utot'*Btot);
            Mdown = Rdown\eye(6*P);

            U2 = Mup*RandVel;
            U1 = Mdown*RandVel;
        end

        RFD = (U2-U1)/delta;

        %RFD = zeros(6*P,1);

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
        %% Solve with dense sqrt of M. The naive and dense way of solving the problem.
        S = stokes_SLP_mat(rvec_in,rvec_out);
        if solve_dense
            %U = solve_mobility(q,rvec_in, rvec_out,Fvec,opt);
            [UT, ~, ~] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,Tblock,Fvec,0, opt);
    
            %build sqrt of M densely:
            [Ytot,Utot] = getPseudoFactors(S,1e-13,0);
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

            %% Solve instead with blockdiag(T)*S iteratively
            A = -T*S;
            %A = -Tdiag*S; %the weights will cancel out
            Asym = M/(4*pi)*(A+A')/2;
            %sqrtA = chol(Asym); %results in very large errors!
            
           %  [Vs,Ds] = eig(Asym);
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
            dW = randn(3*N*P,1); %different noise here for different N!
            % Udense = (R3\eye(P*6))*(Fvec+Ktot'*Ainv*sqrtA*sqrt(2/opt.dt)*dW);
            %noise = sqrtA*dW;

            %% Terrible convergence for Lanczos without precond
            %[noise,~,iters] = KrylovSqrtMsing((A+A')/2,dW,tol_krylov);
            %% Set up Lanczos with preconditioning (implemented in slow way)
            B_precond = C*Asym*C';    
            [noise,iter_errv5,iters] = KrylovSqrtMsing(B_precond,dW,tol_krylov,maxit,Cplus);
            iter_hist(i) = iters;
    %% Solve mobility problem
           % [U, iters, lambda_norm] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,Tblock,Fvec, 1,opt, sqrtA,dW);
            [U, iters, lambda_norm] = solve_brownian_mobility(q,rvec_in,rvec_out,Y,UU,LL,Kin,Tblock,Fvec, 1,opt, [],noise);
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



function [rinUp,routUp,transVec] = updateNodes(rin,rout,P,N,M,vel)

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
