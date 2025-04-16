function res = matvec_mobility_refine(tau,q,rin_base_c,rout_base_c,rin_base_f,...
    rout_base_f,rin_c,rout_c,UU_c,Y_c,Lc,UU_f,Y_f,Lf,vars)
%matvec_moility maps unknowns coefficient vector tau at the particle
%surfaces to flow velocities at the surfaces by first mapping to proxy
%source strengths

P = size(q,1); %number of particles

Mc = size(rout_base_c,1); %points per particle on surface, coarse grid
Nc = size(rin_base_c,1); %points per particle on proxy surface, coarse grid 

Mf = size(rout_base_f,1); %points per particle on surface, fine grid
Nf = size(rin_base_f,1); %points per particle on proxy surface, fine grid 
 
%First, map density for all particles
lambda_stokes = zeros(3*P*Nc,1);
lambda_stokes_cut = zeros(3*P*Nc,1);
lambda_stokes_fine = zeros(3*P*Nf,1);

for i = 1:P
    %Retrieve self evaluation blocks and transform density on surface to
    %density on proxy surface
    if vars.ellipsoid 
% 
          warning('Not yet implemented')
%         Ri = R{i}; % Rotation matrix for the particle 
%     
%         step0 = rotate_vector(tau((i-1)*3*M+1:M*i*3),Ri');
%         step1 = UU*step0; 
%         tau_mapped1 = Y*step1;
%         tau_mapped2 = rotate_vector(tau_mapped1,Ri);
% 
%         % Need to project off any contribution to force and torque
%         lambda_stokes(3*(i-1)*N+1:3*i*N) = tau_mapped2-rotate_vector(L*tau_mapped1,Ri);        

    else
        step1 = UU_c*tau((i-1)*3*Mc+1:Mc*i*3);
        tau_mapped = Y_c*step1; %this is the proxy density in proxy sources
        % Need to project off any contribution to force and torque
        lambda_stokes(3*(i-1)*Nc+1:3*i*Nc) = tau_mapped-Lc*tau_mapped; 

        step1 = UU_c*tau((i-1)*3*Mc+P*Mc*3+1:Mc*i*3+P*Mc*3);
        tau_mapped = Y_c*step1; %this is the proxy density in proxy sources
        % Need to project off any contribution to force and torque
        lambda_stokes_cut(3*(i-1)*Nc+1:3*i*Nc) = tau_mapped-Lc*tau_mapped; 

        %Map also the fine density
        step1 = UU_f*tau((i-1)*3*Mf+P*Mc*6+1:Mf*i*3+P*Mc*6);
        tau_mapped = Y_f*step1; %this is the proxy density in proxy sources
        % Need to project off any contribution to force and torque
        lambda_stokes_fine(3*(i-1)*Nf+1:3*i*Nf) = tau_mapped-Lf*tau_mapped; 
        


    end

end

%% Evaluate flow field
%Do one call to FMM (or direct evaluation)
% NB, later, this should include also the evaluation from image sources!
res1 = getFlow(lambda_stokes,rin_c,rout_c,vars);
res2 = getTruncField(P,q,rin_base_c,rout_base_c,lambda_stokes_cut,vars.Lcut);
res3 = getTruncField(P,q,rin_base_f,rout_base_f,lambda_stokes_fine,vars.Lcut);
res = [res1; res2; res3];
res = res+tau;


%Correct self evaluation. Subtract interaction 
for i = 1:P
    %subtract self constribution
    rin_i = rin_c(Nc*(i-1)+1:Nc*i,:);  
    rows_i = (i-1)*Mc+1:i*Mc;
    targ = rout_c(rows_i,:); 

    stoklet = reshape(lambda_stokes(3*(i-1)*Nc+1:3*i*Nc),3,[]);        
    U = SE0P_Stokeslet_direct_full_ext_mex(rin_i, stoklet', struct('eval_ext_x', targ));
    U_stok = U';
    res(3*(i-1)*Mc+1:3*i*Mc) = res(3*(i-1)*Mc+1:3*i*Mc)-1/8/pi.*U_stok(:);

    stoklet = reshape(lambda_stokes_cut(3*(i-1)*Nc+1:3*i*Nc),3,[]);        
    U = SE0P_Stokeslet_direct_full_ext_mex(rin_i, stoklet', struct('eval_ext_x', targ));
    U_stok = U';
    res(3*(i-1)*Mc+1+3*Mc*P:3*i*Mc+3*Mc*P) = res(3*(i-1)*Mc+1+3*Mc*P:3*i*Mc+3*Mc*P)-1/8/pi.*U_stok(:);

    rin_i = rin_base_f+q(i,:);
    targ = rout_base_f+q(i,:);

    stoklet = reshape(lambda_stokes_fine(3*(i-1)*Nf+1:3*i*Nf),3,[]);        
    U = SE0P_Stokeslet_direct_full_ext_mex(rin_i, stoklet', struct('eval_ext_x', targ));
    U_stok = U';
    res(3*(i-1)*Mf+1+6*Mc*P:3*i*Mf+6*Mc*P) = res(3*(i-1)*Mf+1+6*Mc*P:3*i*Mf+6*Mc*P)-1/8/pi.*U_stok(:);


    %% Need to also subtract blocks corresponding to images, this part is to be written
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
%         res(3*(i-1)*M+1:3*i*M) = res(3*(i-1)*M+1:3*i*M)-1/8/pi*U_stok(:)+...
%             -1/8/pi*U_rot(:)-1/4/pi*U_pot(:)+tau((i-1)*M*3+1:M*i*3);
%     else
%             %subtract and add the identity blocks
%         if vars.precond
%             weight_i = weights((i-1)*M*3+1:i*M*3);
%             res(3*(i-1)*M+1:3*i*M) = res(3*(i-1)*M+1:3*i*M)-1/8/pi*weight_i.*U_stok(:); %+...
%            % tau((i-1)*M*3+1:M*i*3);
%         else
%             res(3*(i-1)*M+1:3*i*M) = res(3*(i-1)*M+1:3*i*M)-1/8/pi.*U_stok(:);
%         end
%     end




end