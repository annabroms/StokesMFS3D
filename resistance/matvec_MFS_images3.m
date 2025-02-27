function res = matvec_MFS_images3(tau,rin,rout,rimage,weights,q,Uii,Yii,vars,R)
% construct matvec

N = size(q,1); %number of particles

N_large = size(rout,1)/N; %points per particle on outer grids
N_small = size(rin,1)/N; %points per particle on proxy surface
N_image = size(rimage,1)/N; %Number of images per particle. NB: For now assumed equal on all particles


%alpha = sum(vars.sing_vec(2:end)); %how many types of singularities to deal with 

 
%Precomputation only for a single paricle
if ~vars.image
    UU = Uii{1};
    Y = Yii{1};
end


%First, map density for all particles
tau_stokes = zeros(3*N*N_small,1);
%tau_image = zeros(3*N*N_image*(vars.beta-1),1);
%tau_rot = zeros(3*N*N_image,1);
%tau_pot = zeros(3*N*N_image,1);
if vars.profile
    memorygraph('label','apply precond in matvec');
end

for i = 1:N
    %Retrieve self evaluation blocks

    if vars.split_precond
        %Preconditioning split for the image block and the rest of the
        %source to target matrix
        UU = Uii{i,1};
        Y = Yii{i,1};


    else
        % if vars.image || vars.ellipsoid %precomputation done for every particle
        %     UU = Uii{i};
        %     Y = Yii{i};
        % end
        %R_large = kron(eye(N_large),R{i});
        %R_small = kron(eye(N_small),R{i});
        if vars.ellipsoid 
            Ri = R{i};
        
            step0 = rotate_vector(tau((i-1)*3*N_large+1:N_large*i*3),Ri');
            step1 = UU*step0;
            tau_mapped = rotate_vector(Y*step1,Ri);
        else
            step1 = UU*tau((i-1)*3*N_large+1:N_large*i*3);
            tau_mapped = Y*step1;
        end


        %step1B = UU*(R_large'*tau((i-1)*3*N_large+1:N_large*i*3));
        

        %rotate_vector(x_gmres((i-1)*3*M+1:i*3*M),R{i}')))

        %Apply QR
        %step1 = L'*(R_large'*tau((i-1)*3*N_large+1:N_large*i*3));        
        %tau_mapped = R_small*(U2(:,1:3*N_small)\step1);

        
        %tau_mapped = Y*step1; %this is the mapped density for this particle to throw in to the kernel
        %tau_mapped = R_small*(Y*step1);
        
        tau_stokes(3*(i-1)*N_small+1:3*i*N_small) = tau_mapped(1:3*N_small);

        %is using LU-decomp instead.


        %tau_image(3*(i-1)*(vars.beta-1)*N_image+1:3*i*N_image*(vars.beta-1)) = tau_mapped(3*N_small+1:end);
       % tau_rot(3*(i-1)*N_image+1:3*i*N_image) = tau_mapped(3*N_small+1:3*N_small+3*N_image);
        %tau_pot(3*(i-1)*N_image+1:3*i*N_image) = tau_mapped(3*N_small+3*N_image+1:end);
    end
end

res = tau; 

if vars.profile
    memorygraph('label','compute FMM');
end

%Do one call to FMM (or direct evaluation)
if vars.fmm && ~vars.image
    nd = 1;
    srcinfo.nd = nd;
    
    ifppreg = 0;
    ifppregtarg = 1;
    
    srcinfo.sources = rin';
    srcinfo.stoklet = reshape(tau_stokes,3,[]);
    
    targ = rout';  
    %eps = 1e-6; % was -6
    eps = vars.eps; 
    
    U = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg);    
    %U = st3ddir(srcinfo,targ,ifppregtarg); %Try to use this one
    if vars.precond
        res = res + 1/4/pi*weights.*U.pottarg(:);
    else
        res = res + 1/4/pi*U.pottarg(:);
    end
    
    %Do this with the routine from SE instead. Better? 
    % srcinfo.stoklet = reshape(tau_stokes,3,[]);
    % U = SE0P_Stokeslet_direct_full_ext_mex(rin, srcinfo.stoklet', struct('eval_ext_x', targ'));
    %    % instead for comparison
    % U = U';
    % res = res + 1/8/pi*U(:);
    % U2 = 1/8/pi*U(:); %for debugging only


    clear U srcinfo;
else
    targ = rout; 
    srcinfo.stoklet = reshape(tau_stokes,3,[]);
    U = SE0P_Stokeslet_direct_full_ext_mex(rin, srcinfo.stoklet', struct('eval_ext_x', targ));
       % instead for comparison
    U = U';
    res = res + 1/8/pi*U(:);

    %% Add image contributions
    %% First from rotlets
    if vars.image
    
        srcinfo.rotlet = reshape(tau_rot,3,[]);
    
        U = SE0P_Rotlet_direct_full_ext_mex(rimage, srcinfo.rotlet', struct('eval_ext_x', targ));
        % instead for comparison
        U = U';
        res = res+1/8/pi*U(:);
       
    
        %% Add the contribution from  potential dipoles
        srcinfo.potdip = reshape(tau_pot,3,[]);
    
        U = SE0P_Potdip_direct_full_ext_mex(rimage, srcinfo.potdip', struct('eval_ext_x', targ));
        % instead for comparison
        U = U';
        res = res+1/4/pi*U(:);
    end


end

if vars.profile
    memorygraph('label','subtract self_interaction');
end

%Correct self evaluation. Subtract interaction 
for i = 1:N
    rin_i = rin(N_small*(i-1)+1:N_small*i,:);
    %srcinfo.sources = rin_i';   
    stoklet = reshape(tau_stokes(3*(i-1)*N_small+1:3*i*N_small),3,[]);    
    
    rows_i = (i-1)*N_large+1:i*N_large;
    %targ = rout(rows_i,:)'; 

    targ = rout(rows_i,:);   
    %eps = 1e-6;
    
    %U = stfmm3d(eps,srcinfo,ifppreg,targ,ifppregtarg);    
    %U = st3ddir(srcinfo,targ,ifppregtarg); % Direct evaluation
    %from FMM but much slower than the counterpart from SE.

    U = SE0P_Stokeslet_direct_full_ext_mex(rin_i, stoklet', struct('eval_ext_x', targ));
   % instead for comparison
    U_stok = U';
    %res(rows) = res(rows)+1/4/pi*U.pottarg(:);
    %res(rows) = res(rows)+1/8/pi*U(:);
   
    %U = st3ddir(srcinfo,targ,ifppregtarg); %Small self-evaluation blocks, use direction summation

    %% Need to also subtract blocks corresponding to images
    if vars.image
        %Start with rotlets
        rim_i = rimage(N_image*(i-1)+1:N_image*i,:);
        srcinfo.rotlet = reshape(tau_rot(3*(i-1)*N_image+1:3*i*N_image),3,[]);
        U = SE0P_Rotlet_direct_full_ext_mex(rim_i, srcinfo.rotlet', struct('eval_ext_x', targ));
        % instead for comparison
        U_rot = U';

        %Contribution from potential dipoles
        srcinfo.potdip = reshape(tau_pot(3*(i-1)*N_image+1:3*i*N_image),3,[]);

        U = SE0P_Potdip_direct_full_ext_mex(rim_i, srcinfo.potdip', struct('eval_ext_x', targ));
        % instead for comparison
        U_pot = U';

        res(3*(i-1)*N_large+1:3*i*N_large) = res(3*(i-1)*N_large+1:3*i*N_large)-1/8/pi*U_stok(:)+...
            -1/8/pi*U_rot(:)-1/4/pi*U_pot(:)+tau((i-1)*N_large*3+1:N_large*i*3);
    else
            %subtract and add the identity blocks
        if vars.precond
            weight_i = weights((i-1)*N_large*3+1:i*N_large*3);
            res(3*(i-1)*N_large+1:3*i*N_large) = res(3*(i-1)*N_large+1:3*i*N_large)-1/8/pi*weight_i.*U_stok(:); %+...
           % tau((i-1)*N_large*3+1:N_large*i*3);
        else
            res(3*(i-1)*N_large+1:3*i*N_large) = res(3*(i-1)*N_large+1:3*i*N_large)-1/8/pi.*U_stok(:);
        end
    end








%     %add images if these are present
%     if N_image
%         rows = setdiff(1:N_large*N*3,(i-1)*N_large*3+1:N_large*i*3); 
% 
%         %Get block from large matrix   
%         Nik_image = [Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i)];
%         for l = 1:alpha-1
%             Nik_image  = [Nik_image Ker_image(rows,N_image*N*3*l+N_image*3*(i-1)+1:N_image*N*3*l+N_image*3*i)];
%         end
%     
%         res(rows) = res(rows) + Nik_image*tau_image(3*(i-1)*N_image*(vars.beta-1)+1:3*i*N_image*(vars.beta-1)); 
%     end

end

end