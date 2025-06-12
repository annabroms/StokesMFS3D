function [U,iters,lambda_norm,uerr] = solve_mobility_refined(q,Fvec,fmm,Lcut)

P = size(q,1); 

if nargin < 4
    Rg = 0.68; %proxy radius
    N = 700; % approximate number of proxy sources on every particle
elseif nargin < 5
    N = 700;
end
Nc = 100;
Rp_c = 0.3;
% Nc = 50;
% Rp_c = 0.2;
% Nc = 700;
% Rp_c = 0.68;
Nf = 700;
Rp_f = 0.68;

%initialize a bunch of parameters. Do not change if you don't really want.
opt = init_MFS(Nc);
opt.fmm = fmm; 
opt.maxit = 200; %max number of gmres iterations
gmres_tol = 1e-7;
opt.plot = 0; 

optf = init_MFS(Nf);
optf.plot = 0; 

opt.Lcut = Lcut;
opt.ellipsoid = 0; 

% Inner proxy surface, outer collocation grid for single sphere using spherical design nodes
[rin_c,rout_c,~,~] =  getGrids(Rp_c,opt); 
[rin_f,rout_f,~,~] =  getGrids(Rp_f,optf); 

%Create pseudoinverse of self-interaction matrix, to be used in one-body
%preconditioning
[Y_c,UU_c,Lc,Kin_c,~] = getSL(rin_c,rout_c);
[Y_f,UU_f,Lf,Kin_f,~] = getSL(rin_f,rout_f);

Nc = size(rin_c,1); %number of sources per particle, coarse grid
Mc = size(rout_c,1); %numer of collocation points per particle, coarse grid
Nf = size(rin_f,1); %number of sources per particle, fine grid
Mf = size(rout_f,1); %numer of collocation points per particle, fine grid

rvec_in_c = [];
rvec_out_c = [];
lambda_vec = []; 

%Create grid on every particle. Also create completion source on every
%particle, given force and torque.

for k = 1:P

    %Create grid 
    rvec_in_c = [rvec_in_c; rin_c+q(k,:)];
    rvec_out_c = [rvec_out_c; rout_c+q(k,:)];
    
    %Create right hand side, given forces and torques on the particles
    F = Fvec(6*(k-1)+1:6*(k-1)+3);
    T = Fvec(6*(k-1)+4:6*k);
    lambda_k = getLambda0(F,T,Kin_c);
    lambda_vec = [lambda_vec; lambda_k];
 
end

%Get flow field due to completion source.
uvec = getFlow(lambda_vec,rvec_in_c,rvec_out_c,opt); 
uvec = -uvec;

uvec = [uvec; zeros(size(uvec)); zeros(size(rout_f,1)*P*3,1)]; %extended system now, as we have the dummy variables too


% figure()
% scatter3(rvec_in_c(:,1),rvec_in_c(:,2),rvec_in_c(:,3),'b')
% hold on
% scatter3(rvec_out_c(:,1),rvec_out_c(:,2),rvec_out_c(:,3),'k')

%% Solve for source strengths
% xx = zeros(size(uvec)); 
% for k = 1:size(uvec,1)
%     k
%     xx(:) = 0;
%     xx(k) = 1;
% 
%     cc(:,k) = matvec_mobility_refine(xx,q,rin_c,rout_c,rin_f,...
%         rout_f,rvec_in_c,rvec_out_c,UU_c,Y_c,Lc,UU_f,Y_f,Lf,opt);
% end

% xx = rand(size(uvec));
% test = matvec_mobility_refine(xx,q,rin_c,rout_c,rin_f,...
%          rout_f,rvec_in_c,rvec_out_c,UU_c,Y_c,Lc,UU_f,Y_f,Lf,opt);

[x_gmres,iters,resvec,real_res] = helsing_gmres(@(x) matvec_mobility_refine(x,q,rin_c,rout_c,rin_f,...
    rout_f,rvec_in_c,rvec_out_c,UU_c,Y_c,Lc,UU_f,Y_f,Lf,opt),uvec,size(uvec,1),opt.maxit,gmres_tol);

%chcek residual
abs_res = norm(matvec_mobility_refine(x_gmres,q,rin_c,rout_c,rin_f,...
    rout_f,rvec_in_c,rvec_out_c,UU_c,Y_c,Lc,UU_f,Y_f,Lf,opt)-uvec);


%% Determine velocities and Map back to the sought density in source points
U = zeros(6*P,1);

for i = 1:P
    lambda_i = Y_c*(UU_c*x_gmres((i-1)*3*Mc+1:i*3*Mc));
    lambda_gmres((i-1)*3*Nc+1:i*3*Nc) = lambda_i;
    U(6*(i-1)+1:6*i) = -Kin_c'*lambda_i; 

    lambda_i = Y_c*(UU_c*x_gmres((i-1)*3*Mc+3*Mc*P+1:i*3*Mc+3*Mc*P));
    lambda_gmres((i-1)*3*Nc+P*3*Nc+1:i*3*Nc+P*3*Nc) = lambda_i;
    U(6*(i-1)+1:6*i) = U(6*(i-1)+1:6*i)+Kin_c'*lambda_i; 

    lambda_i = Y_f*(UU_f*x_gmres((i-1)*3*Mf+6*Mc*P+1:i*3*Mf+6*Mc*P));
    lambda_gmres((i-1)*3*Nf+6*Nc*P+1:i*3*Nf+P*6*Nc) = lambda_i;
    U(6*(i-1)+1:6*i) = U(6*(i-1)+1:6*i)-Kin_f'*lambda_i;
end
lambda_norm = norm(lambda_gmres,inf);


%% Check residual
%Get check points, assuming we work with spheres!
test_vel = 0; 

if test_vel
    b = ellipsoid_param(1,1,1);   % baseline object at the origin, aligned
    b = setupsurfquad(b,[46,55]);
    
    rcheck = []; 
    for k = 1:size(q,1)
        x = q(k,:)' + b.x;    % rot then transl, b just for vis
        rcheck = [rcheck; x'];    
    end
    n_check = size(b.x,2);
    
    
    %plot3(rcheck(:,1),rcheck(:,2),rcheck(:,3),'m.');
    
    %Assign velocities at checkpoints
    ucheck = zeros(n_check*3*size(q,1),1); 
    for k = 1:size(q,1)
        Kcheck = getKmat(rcheck(n_check*(k-1)+1:k*n_check,:),q(k,:));
        ucheck((k-1)*3*n_check+1:3*k*n_check) = Kcheck*U((k-1)*6+1:k*6);
        
    end
    
    for i =1:size(q,1)
        densityK_particle = (eye(3*N)-LL)*lambda_gmres(3*(i-1)*N+1:i*3*N)'+lambda_vec(3*(i-1)*N+1:i*3*N);
        densityK(3*(i-1)*N+1:i*3*N) = densityK_particle;
    end
    
    %get flow and comare
    
    ubdry = getFlow(densityK,rvec_in,rcheck,opt);
    uerr_vec = vecnorm(reshape(ucheck-ubdry,3,[]),2,1)/max(vecnorm(reshape(ucheck,3,[]),2,1));
    uerr = max(uerr_vec);
else

    uerr = [];
end



   

end