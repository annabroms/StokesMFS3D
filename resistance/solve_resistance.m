function [Fvec,iters,lambda_norm,err_res] = solve_resistance(q,U,fmm,Rg,N)
%SOVLE_RESISTANCE(q,Uvec,fmm,Rg,N) solves a Stokes resistance problem for
%spherical particles centered at q moving with translational and angular
%velocities in  U (of size 6P, for each particle a translational velocity
%of size 3 x 1 followed by an angular velocity of size 3 x 1). q is of size
%P x 3. The flag fmm sets if FMM3D is in use. Optional params Rg and N set
%radius of proxy surface and number of source points.
% 
%Returns a vector of forces and torques of size 6P with the same format as
%the input velocity vector. iters is the number of iterations required by
%GMRES and lambda_norm the maximum magnitude of the source vector (bad if
%far too large (loss of accuracy). err_res is the estimated max relative residual 
% at the surface.

%% Some parameters needed for the solve
if nargin < 4
    Rg = 0.7;
    Rg = 0.68;
    N = 700;
elseif nargin < 5
    N = 700;
end

opt = init_MFS(N);
opt.Rg = Rg;

opt.fmm = fmm; 
opt.maxit = 600; %max number of GMRES iterations
gmres_tol = 1e-7;

P = size(q,1); 


% Inner proxy surface, outer collocation grid for single sphere using spherical design nodes
[rin,rout,~,~] =  getGrids(opt.Rg,opt);

N = size(rin,1);
M = size(rout,1); 

rvec_in = [];
rvec_out = [];


for k = 1:P
    %Create grid 
    rvec_in = [rvec_in; rin+q(k,:)];
    rvec_out = [rvec_out; rout+q(k,:)];
    pairs(k,2) = size(rout,1);
     
end

%Visuals
visualise = 0; 
if visualise

    [X,Y,Z] = sphere(12);
    
    figure()
    
    r = 1
    for k = 1:size(q,1); 
        %surf(r*Xs+q(k,1),r*Ys+q(k,2),r*Zs+q(k,3),'EdgeColor','flat');
        hSurface=surf(r*X+q(k,1),r*Y+q(k,2),r*Z+q(k,3));
        set(hSurface,'FaceColor',[1 0 0],'FaceAlpha',0.9,'FaceLighting','gouraud','EdgeColor','none')
        
        %surf(r*Xs+q(k,1),r*Ys+q(k,2),r*Zs+q(k,3),'EdgeColor','blue','FaceColor','blue');
        hold on
        axis equal
    end
    camlight
end


%% Assign RHS in resistance problem
Kout = getKmat(rout,[0,0,0]);
%Maybe a stuct with all the options are needed here
for k = 1:P
    u_bndry((k-1)*3*M+1:3*k*M) = Kout*U((k-1)*6+1:k*6);
end


%% compute preconditioning. It's enough to do this for a single particle 
%as it's assumed that everyone has the same shape.

[Yk,UUk] = precond_ellipsoid_resistance(rin,rout);

%The format is used to prepare for the case when different shapes are is
%use
Y{1} = Yk;
UU{1} = UUk; 
disp('Done preconditioning')

if opt.profile
    memorygraph('label','start matvec resistance');
end

%Solve problem
[mu_gmres,iters,resvec,real_res] = helsing_gmres(@(x) matvec_MFS_res(x,rvec_in,rvec_out,[],[],q,UU,Y,opt,[]),u_bndry',3*size(rvec_out,1),opt.maxit,gmres_tol,0);

if opt.profile
    memorygraph('label','done matvec resistance, remap and determine force')
end

%% Determine source strengths on proxy sources from the solution at the boundary. 
% Also, determine forces and torques on particles

Fvec = zeros(6*size(q,1),1);
Kin = getKmat(rin,[0,0,0]);
for i = 1:P
    lambda_i = Y{1}*(UU{1}*(mu_gmres((i-1)*3*M+1:i*3*M)));
    lambda_gmres(3*(i-1)*N+1:i*3*N) = lambda_i;
    Fvec(6*(i-1)+1:6*i) = Kin'*lambda_i; 
end


%% Evaluate velocity at new nodes

if opt.profile
    memorygraph('label','done solving resistance, compute velocities')
end

% First, get check points on the surface to determine residuals in 
b = ellipsoid(1,1,1);   % baseline object at the oridin, aligned
b = setupsurfquad(b,[46,55]);

rcheck = []; 
for k = 1:P
    x = q(k,:)' + b.x;    % rot then transl, b just for vis
    rcheck = [rcheck; x']; 
end
n_check = size(b.x,2);
%plot3(rcheck(:,1),rcheck(:,2),rcheck(:,3),'m.');

%Assign velocities at checkpoints (read off bc)
ucheck = zeros(n_check*3*P,1); 
for k = 1:P
    Kcheck = getKmat(rcheck(n_check*(k-1)+1:k*n_check,:),q(k,:));
    ucheck((k-1)*3*n_check+1:3*k*n_check) = Kcheck*U((k-1)*6+1:k*6);    
end

%Check flow by evaluating representation and compare
ubdry = getFlow(lambda_gmres,rvec_in,rcheck,opt);
uerr_vec = vecnorm(reshape(ucheck-ubdry,3,[]),2,1)/max(vecnorm(reshape(ucheck,3,[]),2,1));

err_res = max(uerr_vec); 


lambda_norm = norm(lambda_gmres,inf);



end