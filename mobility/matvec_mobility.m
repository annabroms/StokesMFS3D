function res = matvec_mobility(tau,rin,rout,q,UU,Y,L,vars,R)
%matvec_moility maps unknowns coefficient vector tau at the particle
%surfaces to flow velocities at the surfaces by first mapping to proxy
%source strengths

P = size(q,1); %number of particles

M = size(rout,1)/P; %points per particle on surface
N = size(rin,1)/P; %points per particle on proxy surface
 
%First, map density for all particles
lambda_stokes = zeros(3*P*N,1);

for i = 1:P
    %Retrieve self evaluation blocks and transform density on surface to
    %density on proxy surface
    if vars.ellipsoid 

        Ri = R{i}; % Rotation matrix for the particle 
    
        step0 = rotate_vector(tau((i-1)*3*M+1:M*i*3),Ri');
        step1 = UU*step0; 
        tau_mapped1 = Y*step1;
        tau_mapped2 = rotate_vector(tau_mapped1,Ri);

        % Need to project off any contribution to force and torque
        lambda_stokes(3*(i-1)*N+1:3*i*N) = tau_mapped2-rotate_vector(L*tau_mapped1,Ri);        

    else
        step1 = UU*tau((i-1)*3*M+1:M*i*3);
        tau_mapped = Y*step1; %this is the proxy density in proxy sources
        % Need to project off any contribution to force and torque
        lambda_stokes(3*(i-1)*N+1:3*i*N) = tau_mapped-L*tau_mapped; 
    end

end

%% Evaluate flow field
%Do one call to FMM (or direct evaluation)
% NB, later, this should include also the evaluation from image sources!
res = getFlow(lambda_stokes,rin,rout,vars);
res = res+tau;


%Correct self evaluation. Subtract interaction 
for i = 1:P
    rin_i = rin(N*(i-1)+1:N*i,:);  
    stoklet = reshape(lambda_stokes(3*(i-1)*N+1:3*i*N),3,[]);    
    
    rows_i = (i-1)*M+1:i*M;

    targ = rout(rows_i,:);   

    U = SE0P_Stokeslet_direct_full_ext_mex(rin_i, stoklet', struct('eval_ext_x', targ));
    U_stok = U';
    res(3*(i-1)*M+1:3*i*M) = res(3*(i-1)*M+1:3*i*M)-1/8/pi.*U_stok(:);


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