function [Fvec,it,lambda_norm,rhs_norm,abs_res] = solve_resistance(q,U,fmm,Rg,N)

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
opt.gmres_tol = 1e-7;


% Inner proxy surface, outer collocation grid for single sphere using spherical design nodes
[rin,rout,~,weights] =  getGrids(opt.Rg,opt);

N = size(rin,1);
M = size(rout,1); 

rvec_in = [];
rvec_out = [];


for k = 1:size(q,1)
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


Kout = getKmat(rout,[0,0,0]);


%Maybe a stuct with all the options are needed here
for k = 1:size(q,1)
    u_bndry((k-1)*3*M+1:3*k*M) = Kout*U((k-1)*6+1:k*6);
end
rcheck = [100 100 100]; %some point far away -- necessary at the moment to activate fmm
rplane = [100 100 100]; 





[lambda_gmres,utest,it,res,uplane,cc,resvec,rhs_norm,abs_res] = evaluate_Stokes_BVP_fast(rvec_in,rvec_out,[],[],[],q,u_bndry,rcheck,rcheck,rplane,pairs,opt);


% figure(3)
% semilogy(resvec)
% hold on
% title('Convergence resistance')
% grid on
% 
% rcheck = get_sphdesign(1200)+q(1,:);
% u_total = getVelocityField(rvec_in,q,[],[],rcheck,lambda_gmres,size(q,1),N,0,opt,0);
% Kcheck = getKmat(rcheck+q(1,:),[0,0,0]+q(1,:));
% ucheck = Kcheck*U(1:6);
% ucheck = reshape(ucheck,3,size(rcheck,1))';
% utest_vec = reshape(u_total,3,size(rcheck,1))';
% uerr_velocity = max(vecnorm(utest_vec-ucheck,2,2))./max(vecnorm(ucheck,2,2))



%% Determine forces and torques on every particle
Kin = getKmat(rin,[0,0,0]);
Fvec = zeros(6*size(q,1),1);
for i = 1:size(q,1)    
    lambda_i = lambda_gmres((i-1)*3*N+1:3*i*N); 
    Fvec(6*(i-1)+1:6*i) = Kin'*lambda_i; 
end

% figure(1)
% %clf;
% plot(lambda_gmres)
% hold on
lambda_norm = norm(lambda_gmres,inf);



end