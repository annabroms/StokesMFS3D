function [lambda_gmres,u_gmres,it,res,u_plane,traction1,traction2,cc,resvec,rhs_norm,abs_res] = evaluate_Stokes_BVP_fast(rin,rout,rimage,nimage,weights,q,umat,rcheck,rcheck2,rplane,pairs,vars)
%EVALUTATE_STOKES_BVSP_FAST(rin,rout,q,radius,Q,vars) computes the resistance
%matrix with MFS using the single layer operator (the Stokeslet): rin a
%discretisation of the inner surface, rout a discretisation of the outer.
%Particle center coordinates are given in q, and struct with options vars: 


N = size(q,1); %number of particles
N_small = size(rin,1)/N; %vars.fib_n;
N_image = size(rimage,1)/N;

uvec = umat';
rhs_norm = norm(uvec);
%[lambda_gmres,it,res,cc] = solve_Stokes_BVP_split(rin,rout,rimage,nimage,weights,q,[],uvec(:),vars);
[lambda_gmres,it,res,cc,resvec,abs_res] = solve_Stokes_BVP_fast(rin,rout,rimage,nimage,weights,q,pairs,uvec(:),vars);

if vars.plot
    semilogy(abs(lambda_gmres))
end



disp('Determine velocities')

if vars.profile
    memorygraph('label','determine velocity');
end

% We do this for two different sets of test points 
u_gmres = getVelocityField(rin,rout,q,rimage,nimage,rcheck,lambda_gmres,N,pairs,vars,1);
u_plane = getVelocityField(rin,rout,q,rimage,nimage,rplane,lambda_gmres,N,pairs,vars,0);

if vars.get_traction
    traction1 = getTraction(rin, q, N,rcheck, lambda_gmres);
    traction2 = getTraction(rin, q, N,rcheck2, lambda_gmres); %alternative grid
else
    traction1 = [];
    traction2 = [];
end



if vars.profile
    memorygraph('label','all done');
end

end

function NN = remap_source_target(rin,rimage,nimage,singvec,rcheck,N,N_small,N_image,alpha)
    Nio = generate_stokes_matrix_vec_inout(rin,rcheck);
    Ker_image = getImageKernels(singvec,rimage,nimage,rcheck);

    

    NN = [];
    for i = 1:N
        N_proxy = Nio(:,N_small*3*(i-1)+1:N_small*3*i);
        %N_im = Ker_image(:,N_image*3*(i-1)+1:N_image*3*i);
        if sum(singvec(2:end))>0
            N_im = [Ker_image(:,N_image*3*(i-1)+1:N_image*3*i)];
            for k = 1:alpha-1
                N_im = [N_im Ker_image(:,N_image*N*3*k+N_image*3*(i-1)+1:N_image*N*3*k+N_image*3*i)];
            end
        else
            N_im = [];
        end
    
        NN = [NN N_proxy N_im];
    end
end



function u_total = getVelocityField(rin,rout,q,rimage,nimage,rtest,lambda_gmres,N,pairs,vars,flag)

n_proxy = (size(rin,1)-sum(pairs(:,1)))/N; 

if vars.fmm
    % Use the FMM for evaluation
    nd = 1;
    srcinfo.nd = nd;
    ifppreg = 0;
    ifppregtarg = 1;

    srcinfo.sources = rin';   
    srcinfo.stoklet = reshape(lambda_gmres,3,[]);
    
    targ = rtest';  
    eps = 1e-6; % was -6
    U = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg);    
    u_total = 1/4/pi*U.pottarg(:);

else%if vars.image

    %For debugging
    N = size(q,1); %number of particles
    N_small = size(rin,1)/N; %vars.fib_n;
    N_image = size(rimage,1)/N;


    %NN = remap_source_target(rin,rimage,nimage,singvec,rcheck,N,N_small,N_image,alpha);
    %u_gmres = NN*lambda_gmres;
    %u_precond  = NN*lambda_precond;
    
    % lambda_proxy = zeros(3*N*N_small,1);    
    % lambda_image = zeros(3*N*2*N_image,1);

    %NEW
    lambda_proxy = zeros(3*N*n_proxy,1);
    lambda_image = zeros(3*2*sum(pairs(:,1)),1);

    image_ind = 0;
    proxy_ind = 0;
    proxy_ind2 = 0; 

    for k = 1:N

        N_small = proxy_ind +1:proxy_ind+3*n_proxy+3*pairs(k,1);
        N_image = image_ind+1:image_ind+6*pairs(k,1);
        % 
        lambda_proxy(N_small) = lambda_gmres(proxy_ind2+1:proxy_ind2+3*n_proxy+3*pairs(k,1));
        lambda_image(N_image) = lambda_gmres(proxy_ind2+3*n_proxy+3*pairs(k,1)+1:proxy_ind2+3*n_proxy+9*pairs(k,1));

        %lambda_proxy(3*(k-1)*N_small+1:3*k*N_small) = lambda_gmres(3*(k-1)*(N_small+2*N_image)+1:3*(k-1)*(N_small+2*N_image)+3*N_small);
        %lambda_image(2*3*(k-1)*N_image+1:2*3*k*N_image) = lambda_gmres(3*(k-1)*(N_small+2*N_image)+3*N_small+1:3*k*(N_small+2*N_image));
        if vars.mapped
            waring('Does not work for unequal number of sources and target on different particles')
            index1 = (k-1)*N_small+1:(k-1)*N_small+(N_small-N_image);
            dists(index1) = vars.beta*(1-vecnorm(rin(index1,:)-q(k,:),2,2));
            
            index2 = (k-1)*N_small+(N_small-N_image)+1:k*N_small;
            dists(index2) = vars.alpha*(1-vecnorm(rin(index2,:)-q(k,:),2,2));
            %Redo without scaling for stokeslets
            % dists = ones(size(dists)); 
            % dists(index2) = vars.alpha;
        else
            dists = [];
        end
        proxy_ind = proxy_ind+3*n_proxy+3*pairs(k,1);
        proxy_ind2 = proxy_ind2+3*n_proxy+9*pairs(k,1);
        image_ind = image_ind+6*pairs(k,1);
    end

    %LL = getProjection(rtest,q(1,:));
    % if flag
    %     LL = getProjection(rtest,q);
    %    % LL = getProjection(rout,q);
    % end
    
    [u_proxy,u_image] = proxyImageVelocity(lambda_proxy,lambda_image,dists,rin,q,rimage,nimage,rtest,N,pairs,vars);
    
    if N_image
        u_total = (u_proxy+u_image); 
    else 
        u_total = u_proxy;
    end
   
    vars.gamma = 0; 
    if vars.image && vars.image_nbr
        %u_gmres = ((1+vars.gamma)*eye(size(rcheck,1)*3)-vars.gamma*LL)*(u_proxy+u_image); 
               

        if flag && vars.gamma
            N_large = size(rtest,1)/N;
            extra = zeros(size(u_total)); 
        
            for i = 1:N
                lambda_proxy1 = zeros(N_small*N*3,1); 
                lambda_proxy1((i-1)*N_small*3+1:i*N_small*3) = lambda_proxy((i-1)*N_small*3+1:i*N_small*3);

                lambda_image1 = zeros(2*N_image*N*3,1); %two types of images
                lambda_image1(2*(i-1)*N_image*3+1:2*i*N_image*3) = lambda_image(2*(i-1)*N_image*3+1:2*i*N_image*3);
                
    
                [u_proxy1,u_image1] = proxyImageVelocity(lambda_proxy1,lambda_image1,dists,rin,q,rimage,nimage,rtest,N,pairs,vars);
                Li = LL((i-1)*3*N_large+1:i*3*N_large,(i-1)*3*N_large+1:i*3*N_large);
                Li_size = size(Li); 

                %% ONLY WORKS FOR TWO PARTICLES
                if i == 1
                    Li = [Li zeros(Li_size); zeros(Li_size),eye(Li_size)];
                else
                    Li = [eye(Li_size) zeros(Li_size); zeros(Li_size) Li];
                end
                %A = (eye(size(rtest,1)*3)-vars.gamma*Li)*(u_proxy1+u_image1);
                %S = u_proxy1+u_image1;
                B = vars.gamma*(eye(size(rtest,1)*3)-Li);
                extra = extra + B*(u_proxy1+u_image1);
                u_total = u_total+ B*(u_proxy1+u_image1);
            end
        end

        %u_total = (u_proxy+u_image);
    end

end

end

function T = getTraction(rin, q, N,rcheck, lambda_gmres)
    ifppreg = 0;      % no eval at sources
    
     
    srcinfo.sources = rin';
    srcinfo.stoklet = reshape(lambda_gmres,3,[]);  % stokeslet strength vectors
    
    targ = rcheck';  

    ifppregtarg = 3;  % request everything

    M = size(rcheck,1)./N;

    for k = 1:N
        targnor(:,M*(k-1)+1:M*k) = -(rcheck(M*(k-1)+1:M*k,:)-q(k,:))';
    end

%     figure()
%     quiver3(rcheck(:,1),rcheck(:,2),rcheck(:,3),targnor(1,:)',targnor(2,:)',...
%         targnor(3,:)')
    
    tic;
    U = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg);
    toc;
    nt = M*N;
%    fprintf("%d to %d points with u,p,gradu: done in %.3g s (%.3g tot pts/sec)\n",ns,nt,t,(ns+nt)/t)
    p = U.pretarg;         % pressure (1*ntarg)
    gradu = U.gradtarg;    % grad vel (3*3*ntarg)
    shearstress = gradu + permute(gradu,[2 1 3]);        % gradu + gradu^T
    % stress sigma = -pI + mu*(gradu+gradu^T),   then T = sigma.n    We assume mu=1.
    T = -(ones(3,1)*p).*targnor + squeeze(sum(shearstress .* reshape(kron(ones(3,1),targnor),3,3,nt),1));   % compute all tractions (painful tensor contraction in matlab)
%     i = randi(nt);   % which targ to check
%     %fprintf("targ(i=%d): p=%.15g,   gradu = \n",i,p(i)); disp(gradu(:,:,i))
%     fprintf("targ(i=%d) traction T = \n",i);
%     disp(T(:,i))
    T = T./(4*pi);


    sum(T')*4*pi/M./sum(srcinfo.stoklet')
    K = getKmat(rcheck,q);
    F_trac = K'*T(:)*4*pi/M;
    Kin = getKmat(rin,q);
    F_dens = Kin'*lambda_gmres;

end

function [u_proxy,u_image] = proxyImageVelocity(lambda_proxy,lambda_image,dists,rin,q,rimage,nimage,rtest,N,pairs,vars)

N = size(q,1); %number of particles
%N_small = size(rin,1)/N; %vars.fib_n;
%N_image = size(rimage,1)/N;

singvec = vars.sing_vec;
alpha = sum(vars.sing_vec(2:end));

if vars.mapped
    dists = repmat(dists,3,1);
   % dists = dists';
    f = dists(:);
    %f = ones(size(f));
    stoklet = reshape(f.*lambda_proxy,3,[]);
else

    stoklet = reshape(lambda_proxy,3,[]);
end
u_proxy = 1/8/pi*SE0P_Stokeslet_direct_full_ext_mex(rin, stoklet', struct('eval_ext_x', rtest));
u_proxy = u_proxy';
u_proxy = u_proxy(:);

%should write fast method also for this part
Ker_image = getImageKernels(singvec,rimage,nimage,rtest);
NN = [];

image_ind  = 0; 
for i = 1:N
    %N_proxy = Nio(:,N_small*3*(i-1)+1:N_small*3*i);
    %N_im = Ker_image(:,N_image*3*(i-1)+1:N_image*3*i);

    %Order of the unknowns here? 
    if sum(singvec(2:end))>0
        N_image = image_ind+1:image_ind+3*pairs(i,1);
        %N_im = [Ker_image(:,N_image*3*(i-1)+1:N_image*3*i)];
        N_im = [Ker_image(:,N_image)]; %NEW
        for k = 1:alpha-1
            %N_im = [N_im Ker_image(:,N_image*N*3*k+N_image*3*(i-1)+1:N_image*N*3*k+N_image*3*i)];
            N_im = [N_im Ker_image(:,3*sum(pairs(:,1))+3*sum(pairs(1:i-1,1))+1:3*sum(pairs(:,1))+3*sum(pairs(1:i,1)))];
            %%NEW
        end

        if vars.mapped
            warning('NOT implemented for the case with all particles not using the same discretisation')
            rimage_i = rimage(N_image*(i-1)+1:N_image*i,:);
            if size(rimage_i,1)
                dists = 1-vecnorm(rimage_i-q(i,:),2,2);
            else
                dists = []; 
            end
            f = dists;
           % D = repmat(f.^vars.mapped,1,3);
            D = repmat(f,1,3);
            D = D';
            D = D(:); 
            %Z = ones(size(D(:))); %only hits the potential dipoles
            %DD = diag([Z;D]);
            %DD = diag([D;D]);
            DD = diag([D.^2;D.^3]);
            N_im = vars.alpha*N_im*DD;
        end
    else
        N_im = [];
    end

    NN = [NN N_im];
    image_ind = image_ind + 3*pairs(i,1);
end
u_image = NN*lambda_image; 

end






