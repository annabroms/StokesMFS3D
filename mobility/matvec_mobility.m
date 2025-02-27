function res = matvec_mobility(tau,rin,rout,q,UU,Y,L,vars,R)
% construct matvec for mobility 

N = size(q,1); %number of particles

N_large = size(rout,1)/N; %points per particle on outer grid
N_small = size(rin,1)/N; %points per particle on proxy surface
 
%First, map density for all particles
tau_stokes = zeros(3*N*N_small,1);

for i = 1:N
    %Retrieve self evaluation blocks and transform density
    if vars.ellipsoid 
      


        Ri = R{i};
    
        step0 = rotate_vector(tau((i-1)*3*N_large+1:N_large*i*3),Ri');
        step1 = UU*step0; 
        tau_mapped1 = Y*step1;
        tau_mapped2 = rotate_vector(tau_mapped1,Ri);


        %step1 = UU{i}*tau((i-1)*3*N_large+1:N_large*i*3);
        %tau_mapped = Y{i}*step1; %this is the mapped density for this particle to throw in to the kernel
        %tau_stokes(3*(i-1)*N_small+1:3*i*N_small) = (eye(3*N_small)-R_small*L*R_small')*tau_mapped(1:3*N_small);
        %tau_stokes(3*(i-1)*N_small+1:3*i*N_small) = tau_mapped-R_small*(L*(R_small'*tau_mapped));
        %tau_stokes(3*(i-1)*N_small+1:3*i*N_small) = tau_mapped-rotate_vector(L*tau_mapped,Ri);

        tau_stokes(3*(i-1)*N_small+1:3*i*N_small) = tau_mapped2-rotate_vector(L*tau_mapped1,Ri);
        
        %March 12 NOT CORRECT!
        %tau_stokes(3*(i-1)*N_small+1:3*i*N_small) = tau_mapped2-L*tau_mapped2;

    else
        step1 = UU*tau((i-1)*3*N_large+1:N_large*i*3);
        tau_mapped = Y*step1; %this is the mapped density for this particle to throw in to the kernel
        %tau_stokes(3*(i-1)*N_small+1:3*i*N_small) = (eye(3*N_small)-L)*tau_mapped(1:3*N_small);
        tau_stokes(3*(i-1)*N_small+1:3*i*N_small) = tau_mapped-L*tau_mapped; 
    end

end

res = tau; 

%Do one call to FMM (or direct evaluation)
if vars.fmm 
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

    res = res + 1/4/pi*U.pottarg(:);
    
    
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

%     %% Add image contributions
%     %% First from rotlets
%     
%     srcinfo.rotlet = reshape(tau_rot,3,[]);
% 
%     U = SE0P_Rotlet_direct_full_ext_mex(rimage, srcinfo.rotlet', struct('eval_ext_x', targ));
%     % instead for comparison
%     U = U';
%     res = res+1/8/pi*U(:);
%    
% 
%     %% Add the contribution from  potential dipoles
%     srcinfo.potdip = reshape(tau_pot,3,[]);
% 
%     U = SE0P_Potdip_direct_full_ext_mex(rimage, srcinfo.potdip', struct('eval_ext_x', targ));
%     % instead for comparison
%     U = U';
%     res = res+1/4/pi*U(:);


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
    res(3*(i-1)*N_large+1:3*i*N_large) = res(3*(i-1)*N_large+1:3*i*N_large)-1/8/pi.*U_stok(:);
    %res(rows) = res(rows)+1/4/pi*U.pottarg(:);
    %res(rows) = res(rows)+1/8/pi*U(:);
   
    %U = st3ddir(srcinfo,targ,ifppregtarg); %Small self-evaluation blocks, use direction summation

    %% Need to also subtract blocks corresponding to images
%     if vars.image
%         %Start with rotlets
%         rim_i = rimage(N_image*(i-1)+1:N_image*i,:);
%         srcinfo.rotlet = reshape(tau_rot(3*(i-1)*N_image+1:3*i*N_image),3,[]);
%         U = SE0P_Rotlet_direct_full_ext_mex(rim_i, srcinfo.rotlet', struct('eval_ext_x', targ));
%         % instead for comparison
%         U_rot = U';
% 
%         %Contribution from potential dipoles
%         srcinfo.potdip = reshape(tau_pot(3*(i-1)*N_image+1:3*i*N_image),3,[]);
% 
%         U = SE0P_Potdip_direct_full_ext_mex(rim_i, srcinfo.potdip', struct('eval_ext_x', targ));
%         % instead for comparison
%         U_pot = U';
% 
%         res(3*(i-1)*N_large+1:3*i*N_large) = res(3*(i-1)*N_large+1:3*i*N_large)-1/8/pi*U_stok(:)+...
%             -1/8/pi*U_rot(:)-1/4/pi*U_pot(:)+tau((i-1)*N_large*3+1:N_large*i*3);
%     else
%             %subtract and add the identity blocks
%         if vars.precond
%             weight_i = weights((i-1)*N_large*3+1:i*N_large*3);
%             res(3*(i-1)*N_large+1:3*i*N_large) = res(3*(i-1)*N_large+1:3*i*N_large)-1/8/pi*weight_i.*U_stok(:); %+...
%            % tau((i-1)*N_large*3+1:N_large*i*3);
%         else
%             res(3*(i-1)*N_large+1:3*i*N_large) = res(3*(i-1)*N_large+1:3*i*N_large)-1/8/pi.*U_stok(:);
%         end
%     end




end