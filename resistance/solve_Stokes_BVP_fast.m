function [lambda_gmres,it,real_res,cc,resvec,abs_res] = solve_Stokes_BVP_fast(rin,rout,rimage,nimage,weights,q,pairs,uvec,vars)
%solve_BVP_resistance_fast(rin,rout,nin,q,quat,Q,Niit,U,vars) computes a solution to the resistance problem with MFS using the single layer operator (the Stokeslet): rin a
%discretisation of the inner surface, rout a discretisation of the outer.
%Particle center coordinates are given in q, matrix of augumented
%quadrature weights Q and struct with options vars: ls is an upsampling
%factor of the grid, given GL in theta (pole to pole) and trapz in the
%crossection of the sphere.

%warning('This method does not yet converge properly... ')

N = size(q,1); % number of particles
M1 = size(rin,1)/N;
M2 = size(rout,1)/N;
Ntot = size(rin,1);

alpha = sum(vars.sing_vec(2:end)); %how many types of singularities to deal with

Nim = alpha*size(rimage,1)/N;
M3 = M1+Nim;
nout = ceil(vars.ls_glob*(vars.des_n));

n_proxy = (size(rin,1)-sum(pairs(:,1)))/N;



%compute mapped density
%maxit = 100; 
tol = vars.gmres_tol; 
restart = []; 

if vars.precond == 2
    weights = sqrt(weights);
end

%tic
if vars.profile
    memorygraph('label','precompute')
end
[Uii,Yii,Nio,Ker_image] = matvec_MFS_images(rin,rout,rimage,nimage,weights,q,pairs,vars);

%disp('get SVD')
%tic
%[Uii,Yii,Nio,Ker_image,cc] = precond_MFS_unequal(rin,rout,rimage,nimage,weights,q,pairs,vars);
%toc
%toc
%[B2,NiiT2] = matvec_MFS(rin,rout,q,NiiT,vars);
%x = gmres(@(x) matvec_MFS2(x,rin,rout,q,NiiT,vars),uvec,[],tol,maxit);
%x2 = gmres(B,uvec,[],tol,maxit);
% disp('B matrix computed, get conditioning')
% tic
% c = cond(B);
% toc
% tic
% minE = min(abs(eig(B)));
% toc
% disp('Solve with gmres')

cc = nan;


minE = nan; 

%
% A = build_matrix(rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,vars);
% 
% tic 
% cc = skeel(A)
% toc

if vars.precond
    %if solve with preconditioning from the left
%     for i = 1:N
%         WW = repmat(weights((i-1)*M2+1:i*M2),1,3);
%         WW = WW';
%        % uvec(3*(i-1)*M2+1:3*i*M2) = diag(WW(:))*uvec(3*(i-1)*M2+1:3*i*M2);
%         uvec(3*(i-1)*M2+1:3*i*M2) = WW(:).*uvec(3*(i-1)*M2+1:3*i*M2);
%     end

    WW = repmat(weights,1,3);
    WW = WW';
    uvec = WW(:).*uvec;

else
    WW = [];
end

if vars.profile
    memorygraph('label','solve system');
end

%x = ones(size(uvec));
%res = matvec_MFS_images2(x,rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,pairs,vars);


if vars.fmm 
    [x_gmres,it,resvec,real_res] = helsing_gmres_mv(@(x) matvec_MFS_images3(x,rin,rout,rimage,WW(:),q,Uii,Yii,vars),uvec,3*size(rout,1),vars.maxit,tol,1);
    abs_res = norm(matvec_MFS_images3(x_gmres,rin,rout,rimage,WW(:),q,Uii,Yii,vars)-uvec);
else
    [x_gmres,it,resvec,real_res] = helsing_gmres_mv(@(x) matvec_MFS_images_unequal(x,rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,pairs,vars),uvec,3*size(rout,1),vars.maxit,tol,1);
    abs_res = norm(matvec_MFS_images_unequal(x_gmres,rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,pairs,vars)-uvec);

  %  [x_gmres,it,resvec,real_res] = helsing_gmres_mv(@(x) matvec_MFS_images2(x,rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,vars),uvec,3*size(rout,1),vars.maxit,tol);
  %  abs_res = norm(matvec_MFS_images2(x_gmres,rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,vars)-uvec);
end
%[x_gmres2,flag,relres,iter,resvec] = gmres(@(x) matvec_MFS_images2(x,rin,rout,rimage,q,Nio,Ker_image,Uii,Yii,vars),uvec,restart,tol,maxit);

%res = matvec_MFS_images2(x_gmres,rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,vars)-A*x_gmres;
if vars.plot
    figure(77)
    semilogy(resvec)
    grid on 
end

 disp('Gmres done, postprocess')
% x_direct = B\uvec;

% figure(55)
% semilogy(abs(x_direct),'o');
% hold on
% semilogy(abs(x_gmres),'.');


%transform density back
%lambda_gmres = zeros(M3*N*3,1); 
lambda_gmres = zeros(Ntot,1); 
% lambda_direct = zeros(M3*N*3,1); 

% figure(102)
% plot(x);
% hold on
% title('Mapped density','interpreter','latex')
if vars.profile
    memorygraph('label','GMRES done, postprocess')
end

colloc_ind = 0;
source_ind = 0;

tic 
for i = 1:N %


    if vars.image
        if vars.split_precond

            lambda_gmres((i-1)*3*M3+1:(i-1)*3*M3+3*vars.fib_n) = Yii{1,1}*(Uii{1,1}*x_gmres((i-1)*3*M2+1:(i-1)*3*M2+3*nout));
            lambda_gmres((i-1)*3*M3+3*vars.fib_n+1:i*3*M3) = Yii{i,2}*(Uii{i,2}*x_gmres((i-1)*3*M2+3*nout+1:i*3*M2));
        else
            %lambda_gmres((i-1)*3*M3+1:i*3*M3) = Yii{i}*(Uii{i}*x_gmres((i-1)*3*M2+1:i*3*M2));
            %if different number of points on the different particles
            
            N_small = source_ind+1:source_ind+n_proxy*3+9*pairs(i,1);
            N_large = colloc_ind +1:colloc_ind+3*pairs(i,2);
            lambda_gmres(N_small) = Yii{i}*(Uii{i}*x_gmres(N_large));
            colloc_ind = colloc_ind+3*pairs(i,2);
            source_ind = source_ind+3*n_proxy+9*pairs(i,1);
        end
    else
        lambda_gmres((i-1)*3*M3+1:i*3*M3) = Yii{1}*(Uii{1}*x_gmres((i-1)*3*M2+1:i*3*M2));
    end
end
toc
disp('Apply SVD')


%x_orig is the computed density

% figure(56)
% semilogy(abs(lambda_direct),'o');
% hold on
% semilogy(abs(lambda_gmres),'.');



end

function A = build_matrix(rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,vars)
warning('CANNOT TAKE DIFFERENT NUMBER OF SOURCES AND COLLOCATION POINTS PER PARTICLE')
N = size(q,1); %number of particles

N_large = size(rout,1)/N; %points per particle on outer grid
N_small = size(rin,1)/N; %points per particle on proxy surface
N_image = size(rimage,1)/N; %Number of images per particle. NB: For now assumed equal on all particles

alpha = sum(vars.sing_vec(2:end)); %how many types of singularities to deal with
nout = ceil(vars.ls_glob*(vars.fib_n));

if vars.image
    Nik_image = zeros(3*N_large*(N-1),alpha*3*N_image);
end

%Precomputation only for a single paricle
if ~vars.image
    UU = Uii{1};
    Y = Yii{1};
end

A = eye(3*N_large*N,3*N_large*N);

%Create scaling structure
% if vars.mapped
%     for i = 1:N
%         rin_i = rin(N_small*(i-1)+1:N_small*i,:);
%         rimage_i = rimage(N_image*(i-1)+1:N_image*i,:);
%         dists = 1-vecnorm(rimage_i-q(i,:),2,2);
%         f = dists;
% 
%         %D = repmat(f.^vars.mapped,1,3);
%         D = repmat(f,1,3);
%         D_im = D';
%         D_im((i-1)*3*N_image+1:i*N_image) = D(:); 
% 
%         dists = 1-vecnorm(rin_i-q(i,:),2,2);                
%         f = dists;
%         D = repmat(f,1,3);
%         D = D';
%         D(end-3*N_image+1:end) = vars.alpha*D(end-3*N_image+1:end);
%         D_prox((i-1)*3*N_small:i*3*N_small) = D(:); 
%     end
%     %D = ones(size(D));
%     D(end-3*N_image+1:end) = vars.alpha*D(end-3*N_image+1:end);
%     DD = diag(D);
%     N_proxy = N_proxy*DD; 
% end


if vars.precond == 0
    weights = ones(size(weights));
end
   
for i = 1:N
    %Retrieve self evaluation blocks
   
    if vars.image %precomputation done for every particle
        UU = Uii{i};
        Y = Yii{i};
    end
    AiiT = Y*UU;

    
        

%    % Pick all the other indices not on patch
    if vars.split_precond
        UU = Uii{i,2};
        Y = Yii{i,2};
        AiiT = Y*UU;
        rows = setdiff(1:N_large*N*3,(i-1)*N_large*3+3*nout+1:N_large*i*3);
        w_rows = setdiff(1:N_large*N,(i-1)*N_large+nout:N_large*i);

        if N_image && sum(vars.sing_vec(2:end))>0
            %Nik_image(:,1:3*N_image) = [Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i)];
            Nik_image = [Nio(rows,N_small*(i-1)*3+(N_small-N_image)*3+1:N_small*i*3) Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i)];
            %Nik_image(:,1:3*N_image)
            for l = 1:alpha-1
                Nik_image  = [Nik_image Ker_image(rows,N_image*N*3*l+N_image*3*(i-1)+1:N_image*N*3*l+N_image*3*i)];
           
            end
               
            
        end

        A(rows,N_large*3*(i-1)+3*nout+1:N_large*3*i) = Nik_image*AiiT;

        %Now, do the same thing for the proxy points 
        UU = Uii{1,1};
        Y = Yii{1,1};
        AiiT = Y*UU;
        rows = setdiff(1:N_large*N*3,(i-1)*N_large*3+1:(i-1)*N_large*3+3*nout);

        %rows = N_large*(i-1)*3+3*nout+1:N_large*i*3;
        Nii = Nio(rows,N_small*3*(i-1)+1:N_small*3*(i-1)+3*(N_small-N_image));
        A(rows,N_large*3*(i-1)+1:N_large*3*(i-1)+3*nout) = Nii*AiiT; 



    else
%       
        rows = setdiff(1:N_large*N*3,(i-1)*N_large*3+1:N_large*i*3);
        w_rows = setdiff(1:N_large*N,(i-1)*N_large+1:N_large*i);
    
        WW = repmat(weights(w_rows),1,3);
        WW = WW';
    
        %If contribution from image points
        
        % Do this with precomputed matrices: 
        Nik_proxy = Nio(rows,N_small*3*(i-1)+1:N_small*3*i);
    
        if N_image && sum(vars.sing_vec(2:end))>0
            %Nik_image(:,1:3*N_image) = [Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i)];
            Nik_image = Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i);
            %Nik_image(:,1:3*N_image)
            for l = 1:alpha-1
                Nik_image(:,l*3*N_image+1:(l+1)*3*N_image)  = Ker_image(rows,N_image*N*3*l+N_image*3*(i-1)+1:N_image*N*3*l+N_image*3*i);
                %Nik_image  = [Nik_image Ker_image(rows,N_image*N*3*l+N_image*3*(i-1)+1:N_image*N*3*l+N_image*3*i)];
            end
               
            
        end
        
        if vars.mapped
            rin_i = rin(N_small*(i-1)+1:N_small*i,:);
    
            if N_image
                rimage_i = rimage(N_image*(i-1)+1:N_image*i,:);
            
                dists = 1-vecnorm(rimage_i-q(i,:),2,2);
                f = dists;
                
                %D = repmat(f.^vars.mapped,1,3);
                D = repmat(f,1,3);
                D = D';
                D = D(:); 
                DD_im = vars.alpha*diag([D.^2;D.^3]);
            else
                DD_im = [];
            end
        
            dists = 1-vecnorm(rin_i-q(i,:),2,2);                
            f = dists;
            D = repmat(f,1,3);
            D = D';
            D = D(:);
            D(1:3*(N_small-N_image)) = vars.beta*D(1:3*(N_small-N_image));
            D(end-3*N_image+1:end) = vars.alpha*D(end-3*N_image+1:end);
            DD_prox = diag(D); 
        
        
            
            
            A(rows,N_large*3*(i-1)+1:N_large*3*i) = diag(WW(:))*[Nik_proxy*DD_prox Nik_image*DD_im]*AiiT;
        else
            A(rows,N_large*3*(i-1)+1:N_large*3*i) = diag(WW(:))*[Nik_proxy Nik_image]*AiiT;
        end

    end

              
 
end


end

function c = estCondNbr(rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,vars)
N = 1000;
rng(6); 
A = hilb(50);
%A = wilkinson(500);

x = randn(50,1);
%x = randn(3*size(rout,1),1);
r = A*x; 
%r = matvec_MFS_images2(x,rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,vars);
rnorm = norm(r);
xnorm = norm(x); 
for i = 1:N
    dx = 1e-5*randn(3*size(rout,1),1);
    dx = 1e-4*randn(50,1);
    x2 = x+dx;
    dr = A*x2; 
    %dr = matvec_MFS_images2(x2,rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,vars);
    res(i) = (norm(r-dr)/rnorm)/(norm(dx)/xnorm);
    %res(i) = (rnorm/norm(r-dr))/(norm(dx)/norm(x2));
end

c = mean(res); 
end

