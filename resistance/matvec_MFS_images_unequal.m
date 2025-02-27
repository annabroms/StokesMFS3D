function res = matvec_MFS_images_unequal(tau,rin,rout,rimage,weights,q,Nio,Ker_image,Uii,Yii,pairs,vars)
% construct matvec

N = size(q,1); %number of particles
n_proxy = (size(rin,1)-sum(pairs(:,1)))/N;

%N_large = size(rout,1)/N; %points per particle on outer grid
%N_small = size(rin,1)/N; %points per particle on proxy surface
%N_image = size(rimage,1)/N; %Number of images per particle. NB: For now assumed equal on all particles

alpha = sum(vars.sing_vec(2:end)); %how many types of singularities to deal with

%nout = ceil(vars.ls_glob*(vars.fib_n));
 
direct = vars.dense; 

% if vars.image
% 
%     Nik_image = zeros(3*N_large*(N-1),alpha*3*N_image);
% end

%Precomputation only for a single paricle
if ~vars.image
    UU = Uii{1};
    Y = Yii{1};
end

res = tau; %identity mapping

image_ind = 0;
colloc_ind1 = 0;
colloc_ind2 = 0;
proxy_ind = 0;
   
for i = 1:N
    %Retrieve self evaluation blocks
    if vars.split_precond
        %Preconditioning split for the image block and the rest of the
        %source to target matrix

        %For the part that is constant:
        UU = Uii{1,1};
        Y = Yii{1,1};
        
        step1 = UU*tau((i-1)*3*N_large+1:(i-1)*3*N_large+3*nout);
        tau_proxy = Y*step1;

        % From the images
        UU = Uii{i,2};
        Y = Yii{i,2};
        step1 = UU*tau((i-1)*3*N_large+3*nout+1:i*3*N_large);
        tau_image = Y*step1;
       % tau_mapped = [tau_proxy; tau_image];


        %NB: should be modified if we do two different pseudoinverses!
        %Pick rows corresponding to proxy sources and everyone else! 
        rows = setdiff(1:N_large*N*3,(i-1)*N_large*3+3*nout+1:N_large*i*3);
       % rows = (i-1)*N_large*3+1:N_large*(i-1)*3+3*nout; 
      
        
        Nii_image = [Nio(rows,N_small*3*(i-1)+3*(N_small-N_image)+1:N_small*3*i)  Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i)];
        for l = 1:alpha-1
            Nii_image  = [Nii_image Ker_image(rows,N_image*N*3*l+N_image*3*(i-1)+1:N_image*N*3*l+N_image*3*i)];
        end
        res(rows) = res(rows) + Nii_image*tau_image;
    
        %Pick everybody but the proxys on this particle
        rows = setdiff(1:N_large*N*3,(i-1)*N_large*3+1:(i-1)*N_large*3+3*nout);

        %rows = N_large*(i-1)*3+3*nout+1:N_large*i*3;
        Nii = Nio(rows,N_small*3*(i-1)+1:N_small*3*(i-1)+3*(N_small-N_image));
        res(rows) = res(rows) + Nii*tau_proxy;



    else
        if vars.image %precomputation done for every particle
            UU = Uii{i};
            Y = Yii{i};
        end
        N_large = colloc_ind2+1:colloc_ind2+3*pairs(i,2);
        %step1 = UU*tau((i-1)*3*N_large+1:N_large*i*3);
        step1 = UU*tau(N_large);
        tau_mapped = Y*step1; %this is the mapped density for this particle to throw in to the kernel
        %tau_stokes(3*(i-1)*N_small+1:3*i*N_small) = tau_mapped(1:3*N_small);
        %tau_image(3*(i-1)*(vars.beta-1)*N_image+1:3*i*N_image*(vars.beta-1)) = tau_mapped(3*N_small+1:end);
    
    end

% 
%    % Pick all the other indices and multiply through with the pseudoinverse
%    % for particle i
%     
     


    %If contribution from image points
    if vars.split_precond
        %% Add image point contribution
%         rows =  setdiff(1:N_large*N*3,(i-1)*N_large*3+1:(i-1)*N_large*3+3*nout);
%         %w_rows = setdiff(1:N_large*N,(i-1)*N_large+1:(i-1)*N_large+nout);
%         
% %         Nik_image(:,1:6*N_image) = [Nio(rows,N_small*3*(i-1)+3*(N_small-N_image)+1:N_small*3*i)  Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i)];
% %         for l = 1:alpha-1
% %             Nik_image(:,(alpha+1)*3*N_image+1:(alpha+2)*3*N_image)  = [Nik_image Ker_image(rows,N_image*N*3*l+N_image*3*(i-1)+1:N_image*N*3*l+N_image*3*i)];
% %         end
% %         res(rows) = res(rows) + Nik_image*tau_image;
% 
%         %% Add proxy point contribution
%         Nik_proxy = Nio(rows,N_small*3*(i-1)+1:N_small*3*(i-1)+(N_small-N_image)*3);
%         res(rows) = res(rows) + Nik_proxy*tau_proxy;


    else
        %IF ALL HAVE SAME DISCRETISATION:
        %rows = setdiff(1:N_large*N*3,(i-1)*N_large*3+1:N_large*i*3);
        %w_rows = setdiff(1:N_large*N,(i-1)*N_large+1:N_large*i);


        rows = setdiff(1:sum(pairs(:,2))*3,N_large);
        N_rows = colloc_ind1+1:colloc_ind1+pairs(i,2);
        w_rows = setdiff(1:sum(pairs(:,2)),N_rows);


    if direct     
        % Do this with precomputed matrices: SHOULD NOW ONLY WORK IF
        % PARTICLES HAVE THE SAME DISCRETIZATION
        Nik_proxy = Nio(rows,N_small*3*(i-1)+1:N_small*3*i);
        res(rows) = res(rows) + Nik_proxy*tau_mapped(1:3*N_small);

        if N_image && sum(vars.sing_vec(2:end))>0
            %Nik_image(:,1:3*N_image) = [Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i)];
            Nik_image = Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i);
            %Nik_image(:,1:3*N_image)
            for l = 1:alpha-1
                Nik_image(:,l*3*N_image+1:(l+1)*3*N_image)  = Ker_image(rows,N_image*N*3*l+N_image*3*(i-1)+1:N_image*N*3*l+N_image*3*i);
                %Nik_image  = [Nik_image Ker_image(rows,N_image*N*3*l+N_image*3*(i-1)+1:N_image*N*3*l+N_image*3*i)];
            end
        
            res(rows) = res(rows) + Nik_image*tau_mapped(3*N_small+1:end);
        end

    else
       
        targ = rout(w_rows,:); 
        %% Add image contributions
        %% First from rotlets

        if vars.image && vars.image_nbr
            
            %rim_i = rimage(N_image*(i-1)+1:N_image*i,:);
            N_image = image_ind+1:image_ind+pairs(i,1);
            rim_i = rimage(N_image,:);
            
             if vars.mapped
                dists = 1-vecnorm(rim_i-q(i,:),2,2); %how to spread these weights?
                dists = repmat(dists,1,3);
                dists = dists';
                f = dists(:); 
             %end
                % rotlets = f.^vars.mapped.*tau_mapped(3*N_small+1:3*N_small+3*N_image);
                 rotlets = vars.alpha*f.^2.*tau_mapped(3*N_small+1:3*N_small+3*N_image);
            else

                %rotlets = tau_mapped(3*N_small+1:3*N_small+3*N_image);
                rotlets = tau_mapped(3*n_proxy+3*pairs(i,1)+1:3*n_proxy+6*pairs(i,1));
            end

                       
            srcinfo.rotlet = reshape(rotlets,3,[]);
    
            U = SE0P_Rotlet_direct_full_ext_mex(rim_i, srcinfo.rotlet', struct('eval_ext_x', targ));
            % instead for comparison
            U = U';
            %LL = getProjection(targ,q([1:i-1 i+1:N],:));
            
            if vars.precond
                WW = repmat(weights(w_rows),1,3);
                WW = WW';                
               % res(rows) = res(rows)+((1+vars.gamma)*eye(N_large*3)-vars.gamma*LL)*1/8/pi*WW(:).*U(:);
                res(rows) = res(rows)+1/8/pi*WW(:).*U(:);
            else
               % res(rows) = res(rows)+((1+vars.gamma)*eye(N_large*3)-vars.gamma*LL)*1/8/pi*U(:);
                res(rows) = res(rows)+1/8/pi*U(:);
            end
           
    
            %% Add the contribution from  potential dipoles

            if vars.mapped
               % potdip = f.^vars.mapped.*tau_mapped(3*N_small+3*N_image+1:end);
                potdip = vars.alpha*f.^3.*tau_mapped(3*N_small+3*N_image+1:end);
            else
                %potdip = tau_mapped(3*N_small+3*N_image+1:end);
                potdip = tau_mapped(3*n_proxy+6*pairs(i,1)+1:end);
            end


            srcinfo.potdip = reshape(potdip,3,[]);

            
    
            U = SE0P_Potdip_direct_full_ext_mex(rim_i, srcinfo.potdip', struct('eval_ext_x', targ));
            % instead for comparison
            U = U';
            if vars.precond
                %res(rows) = res(rows)+((1+vars.gamma)*eye(N_large*3)-vars.gamma*LL)*1/4/pi*WW(:).*U(:);
                res(rows) = res(rows)+1/4/pi*WW(:).*U(:);
            else
               % res(rows) = res(rows)+((1+vars.gamma)*eye(N_large*3)-vars.gamma*LL)*1/4/pi*U(:);
                res(rows) = res(rows)+1/4/pi*U(:);

            end
        end
       

        %% Add contribution from Stokeslets
        %rin_i = rin(N_small*(i-1)+1:N_small*i,:);  
        N_small = proxy_ind +1:proxy_ind+n_proxy+pairs(i,1);
        rin_i = rin(N_small,:);

        if vars.mapped
               % potdip = f.^vars.mapped.*tau_mapped(3*N_small+3*N_image+1:end);
            dists = 1-vecnorm(rin_i-q(i,:),2,2); %how to spread these weights?
            dists = repmat(dists,1,3);
            dists = dists';
            f = dists(:);
            %f = ones(size(f));
            f(1:3*(N_small-N_image)) = vars.beta*f(1:3*(N_small-N_image));
            f(3*(N_small-N_image)+1:end) = vars.alpha*f(3*(N_small-N_image)+1:end);
            stoklet = f.*tau_mapped(1:3*N_small);
        else
           % stoklet = tau_mapped(1:3*N_small);
            stoklet = tau_mapped(1:3*(n_proxy+pairs(i,1)));
        end


        srcinfo.stoklet = reshape(stoklet,3,[]);
        
        U = SE0P_Stokeslet_direct_full_ext_mex(rin_i, srcinfo.stoklet', struct('eval_ext_x', targ));
       % instead for comparison
        U = U';
        if vars.precond
          %  res(rows) = res(rows)+((1+vars.gamma)*eye(N_large*3)-vars.gamma*LL)*1/8/pi*WW(:).*U(:);
            WW = repmat(weights(w_rows),1,3);
            WW = WW';   
            res(rows) = res(rows)+1/8/pi*WW(:).*U(:);
        else
           % res(rows) = res(rows)+((1+vars.gamma)*eye(N_large*3)-vars.gamma*LL)*1/8/pi*U(:);
            res(rows) = res(rows)+1/8/pi*U(:);
        end

    end              
 
    end

    colloc_ind1 = colloc_ind1+pairs(i,2);
    colloc_ind2 = colloc_ind2+3*pairs(i,2);
    proxy_ind = proxy_ind+n_proxy+pairs(i,1);
    image_ind = image_ind+pairs(i,1);


end