function [U,iters,lambda_norm,uerr] = solve_mobility_refined_RPY(q,Fvec,fmm,Lcut)
%solve with RPY instead. Do shells with spherical designs. Solve without
%preconditioning

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
% Nc = 500;
% Rp_c = 0.5;
Nf = 700;
Rp_f = 0.68;

Rp_c = 0.971774447859209; %for Nc = 100
a_c = 0.109681640780735;
Rp_f = 0.985089125344639; % for Nc = 700
a_f = 0.049918993480785;

% Rp_c = Rp_f;
% Nc = Nf;
% a_c = a_f; 

%initialize a bunch of parameters. Do not change if you don't really want.
opt = init_MFS(Nc);
opt.fmm = fmm; 
opt.maxit = 200; %max number of gmres iterations
gmres_tol = 1e-14;
opt.plot = 0; 

optf = init_MFS(Nf);
optf.plot = 0; 

opt.Lcut = Lcut;
opt.ellipsoid = 0; 

warning('optimization needed')

% Inner proxy surface, outer collocation grid for single sphere using spherical design nodes
[rin_c,rout_c,~,~] =  getGrids(Rp_c,opt); 
[rin_f,rout_f,~,~] =  getGrids(Rp_f,optf); 


Nc = size(rin_c,1); %number of sources per particle, coarse grid
Nf = size(rin_f,1); %number of sources per particle, fine grid


rvec_in_c = [];
rvec_in_f = [];
lambda_fine = []; 
lambda_coarse = []; 

%Create grid on every particle. Also create completion source on every
%particle, given force and torque.

for k = 1:P

    %Create grid 
    rvec_in_c = [rvec_in_c; rin_c+q(k,:)];
    rvec_in_f = [rvec_in_f; rin_c+q(k,:)];
    
 
end

uvec = [zeros(3*Nc*P,1); zeros(3*Nc*P,1); zeros(3*Nf*P,1); -Fvec];

%uvec = [uvec; zeros(size(uvec)); zeros(size(rout_f,1)*P*3,1)]; %extended system now, as we have the dummy variables too


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
%     cc(:,k) = matvec_mobility_RPY(xx,q,rin_c,rout_c,rin_f,...
%         rout_f,rvec_in_c,rvec_out_c,opt);
% end

% xx = rand(size(uvec));
% test = matvec_mobility_refine(xx,q,rin_c,rout_c,rin_f,...
%          rout_f,rvec_in_c,rvec_out_c,UU_c,Y_c,Lc,UU_f,Y_f,Lf,opt);


RPYc = generate_blob_matrix_vec(rvec_in_c,a_c);
RPYcut = getTruncRPY(P,q,rin_c,a_c,Lcut);
RPYfine = getTruncRPY(P,q,rin_f,a_f,Lcut);

[x_gmres,iters,resvec,real_res] = helsing_gmres(@(x) matvec_mobility_refine_RPY(x,q,rin_c,...
    rin_f,RPYc,RPYcut,RPYfine),uvec,size(uvec,1),opt.maxit,gmres_tol);

%chcek residual
abs_res = norm(matvec_mobility_refine_RPY(x_gmres,q,rin_c,...
    rin_f,RPYc,RPYcut,RPYfine)-uvec)


%% Determine velocities 
U = x_gmres(6*Nc*P+3*Nf*P+1:end);
lambda_norm = norm(x_gmres,inf);


%% Check residual
%Get check points, assuming we work with spheres!
test_vel = 1; 

if test_vel
    b = ellipsoid_param(1,1,1);   % baseline object at the origin, aligned
    b = setupsurfquad(b,[46,55]);
    
    rcheck = []; 
    for k = 1:P
        x = q(k,:)' + b.x;    % rot then transl, b just for vis
        rcheck = [rcheck; x'];    
    end
    n_check = size(b.x,2);
    
    
    %plot3(rcheck(:,1),rcheck(:,2),rcheck(:,3),'m.');
    
    %Assign velocities at checkpoints
    ucheck = zeros(n_check*3*size(q,1),1); 
    for k = 1:P
        Kcheck = getKmat(rcheck(n_check*(k-1)+1:k*n_check,:),q(k,:));
        ucheck((k-1)*3*n_check+1:3*k*n_check) = Kcheck*U((k-1)*6+1:k*6);
        
    end
    
    for i =1:P
        densityK_particle = x_gmres(3*(i-1)*Nc+1:i*3*Nc)';
        densityCoarse(3*(i-1)*Nc+1:i*3*Nc) = densityK_particle;

        densityK_particle = x_gmres(3*(i-1)*Nc+1+3*P*Nc:i*3*Nc+3*P*Nc)';
        densityCut(3*(i-1)*Nc+1:i*3*Nc) = densityK_particle;

        densityK_particle = x_gmres(3*(i-1)*Nf+1+6*P*Nc:i*3*Nf+6*P*Nc)';
        densityFine(3*(i-1)*Nf+1:i*3*Nf) = densityK_particle;
    end
    
    %get flow and compare
    
    warning('Ignore the following')
    ucoarse = getFlow(densityCoarse,rvec_in_c,rcheck,opt);
    ucut = getTruncField(P,q,rin_c,b.x',densityCut',Lcut);
    ufine = getTruncField(P,q,rin_f,b.x',densityFine',Lcut);
    ubdry = ucoarse+ucut+ufine;
    ubdry = ucoarse-ucut+ufine;

    uerr_vec = vecnorm(reshape(ucheck-ubdry,3,[]),2,1)/max(vecnorm(reshape(ucheck,3,[]),2,1));
    uerr = max(uerr_vec)
else

    uerr = [];
end



   

end

