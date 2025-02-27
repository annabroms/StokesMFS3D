function [Uii,Yii,SLP,Ker_image] = matvec_MFS_images(rin,rout,rimage,nimage,weights,q,NiiT,vars)
%[B,NiiT,Uii,Yii,Nio,Ker_image] = matvec_MFS_images(rin,rout,rimage,nimage,quat,NiiT,vars)
% construct self-interaction blocks and if wanted, the global source to
% target matrix, divided into the SLP matrix and the matrix for the image points 

N = size(q,1); %number of particles

% For now, assume equal number of source and allocation points for every particle
N_large = size(rout,1)/N; 
N_small = size(rin,1)/N;
N_image = size(rimage,1)/N;
%nout = ceil(vars.ls_glob*vars.fib_n);%+vars.image_nbr*vars.beta*vars.contacts)) %are these extra points important?

alpha = sum(vars.sing_vec(2:end)); %how many types of singularities to deal with
 
%Extract global SLP matrix
if vars.dense
    SLP = generate_stokes_matrix_vec_inout(rin,rout);
    %Extract global Image Kernal matrix
    Ker_image = getImageKernels(vars.sing_vec,rimage,nimage,rout);
else
     SLP = []; 
     Ker_image = []; 
end


% For debugging only
% % [UU,S,V] = svd([SLP Ker_image]);  
% % S = diag(S);
% % figure(54)
% % clf;
% % semilogy(S);
% % c = max(S)/min(S);  %Condition number
% % str = sprintf('Big matrix condition number %1.3e, max sing %1.3e',c,max(S));
% % title(str,'interpreter','latex');

if ~vars.image
    N = 1; %store only preconditioning for a single particle
end
  
%Preparte to store self-interaction blocks
if vars.image && vars.split_precond
    Uii = cell(N,2);
    Yii = cell(N,2);

else
    Uii = cell(N,1);
    Yii = cell(N,1); 
    %Lii = cell(N,1);
end

  % B = eye(3*N_large*N,3*N_large*N);
if vars.image
    last_index = N;
else
    last_index = 1; 
end


    for i = 1:last_index%:N
       
        
        %Take out N for this specific particle
        if vars.dense
            %Extract self evaluation blocks from big global matrix
            rows = N_large*3*(i-1)+1:N_large*3*i;
            N_proxy = SLP(rows,N_small*3*(i-1)+1:N_small*3*i); 
    
            if ~N_image
                N_im = [];
            else
                N_im = [Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i)];
                for k = 1:alpha-1
                    N_im = [N_im Ker_image(rows,N_image*N*3*k+N_image*3*(i-1)+1:N_image*N*3*k+N_image*3*i)];
                end
            end
        else
            %Alternatively, create the blocks here only for particle i
            rin_i = rin(N_small*(i-1)+1:N_small*i,:);
            rout_i = rout(N_large*(i-1)+1:N_large*i,:);
            rimage_i = rimage(N_image*(i-1)+1:N_image*i,:);
            nimage_i = nimage(N_image*(i-1)+1:N_image*i,:);
    
            N_proxy = generate_stokes_matrix_vec_inout(rin_i,rout_i);  
            
            
            N_im = getImageKernels(vars.sing_vec,rimage_i,nimage_i,rout_i);
        end
        if vars.mapped
            rin_i = rin(N_small*(i-1)+1:N_small*i,:);
            rout_i = rout(N_large*(i-1)+1:N_large*i,:);
            rimage_i = rimage(N_image*(i-1)+1:N_image*i,:);
            if size(rimage_i,1)
                dists = 1-vecnorm(rimage_i-q(i,:),2,2);
            else
                dists = [];
            end
            f = 1./dists;
            %D = repmat(f.^vars.mapped,1,3);
            D = repmat(f.^(-1),1,3);
            D = D';
            D = D(:); 
            %Z = ones(size(D(:))); %only hits the potential dipoles
           % DD = diag([Z;D]);
            %DD = diag([D;D]);
            DD = diag([D.^2;D.^3]);
            N_im = vars.alpha*N_im*DD;

            dists = 1-vecnorm(rin_i-q(i,:),2,2);                
            f = dists;
            D = repmat(f,1,3);
            D = D';
            D = D(:); 
            %D = ones(size(D));
            D(end-3*N_image+1:end) = vars.alpha*D(end-3*N_image+1:end);
            D(1:end-3*N_image) = vars.beta*D(1:end-3*N_image);
            DD = diag(D);
            N_proxy = N_proxy*DD; 
        end

        
       
        
        if vars.split_precond
           
            N_p = N_proxy(1:3*nout,1:3*vars.fib_n);
            N_i = [N_proxy(3*nout+1:end,3*vars.fib_n+1:end) N_im(3*nout+1:end,:)];
            
            %Determine the SVD of the two blocks separately
            if i == 1
                [Y,Bi1] = getSVDFactors(N_p,vars.tol(1),vars.plot);
                Uii{i,1} = Bi1';
                Yii{i,1} = Y;
            end
            
            [Y,Bi1]  = getSVDFactors(N_i,vars.tol(2),vars.plot);
            Uii{i,2} = Bi1';
            Yii{i,2} = Y;


        else
            Nii = [N_proxy N_im];
            %Nii = lap3dchgpotmat(rout,rin); %just for debugging --  compare to Alex's code for laplace

           %  body_outer = ellipsoid(1,1,1);   % baseline object at the oridin, aligned
           %  body_outer = setupsurfquad(body_outer,[round(vars.ls*vars.nph),round(vars.ls*vars.nth)]);
           % 
           %  %WW = repmat(body_outer.w,1,3);

           %  WW = repmat(body_outer.w,3,1);
           % % WW = WW';
           %  WW = WW(:);

           %% Test to precondition
           
%            WW = repmat(sqrt(ptArea),1,3);
%            WW = WW';
%            Nii = diag(WW(:))*Nii;

           
           if vars.precond
               WW = repmat(weights((i-1)*N_large+1:i*N_large),1,3);
               WW = WW';
               Nii = diag(WW(:))*Nii;
           end

           %LL = getProjection(rout_i,q(i,:));
            
           %Nii = ((1+vars.gamma)*eye(N_large*3)-vars.gamma*LL)*Nii;

            [Y,Bi1]  = getSVDFactors(Nii,vars.tol,vars.remove_last,vars.plot); 

            Uii{i} = Bi1';
            Yii{i} = Y; 
          %  Lii{i} = LL;
        end
        
        % ONLY IF WE WANT TO EXPLICITLY BUILD B TO INVESTIGATE CONDITIONING
        %THIS IS HOWEVER A BAD IDEA AS B IS NOT STABLY CONSTRUCTED. We then
        %would need to aslo build N globally. 

        %AiiT = Y*(Bi1');
        %Loop through all other blocks in the same column of N and multiply
        %from the right with AiiT
        % ind = [1:i-1 i+1:N];
        % for k = ind
        %     rows = N_large*3*(k-1)+1:N_large*3*k;
        %     Nik_proxy = SLP(rows,N_small*3*(i-1)+1:N_small*3*i);
        %     Nik_image = [Ker_image(rows,N_image*3*(i-1)+1:N_image*3*i)];
        %     for l = 1:alpha-1
        %         Nik_image  = [Nik_image Ker_image(rows,N_image*N*3*l+N_image*3*(i-1)+1:N_image*N*3*l+N_image*3*i)];
        %     end
        % 
        %    % Nik_image = Ker_image(N_large*3*(k-1)+1:N_large*3*k,N_image*3*(i-1)+1:N_image*3*i);
        %     %Nik = [Nik_proxy Nik_image];
        %     %B(rows,N_large*3*(i-1)+1:N_large*3*i) = Nik*AiiT;
        % end       
                    
     
    end

end

function  [Y,Bi1]  = getSVDFactors(N,tol,truncate_last,visualise)

if nargin < 4
    visualise = 0;
end

[UU,S,V] = svd(N);
S = diag(S);
%ra = sum(S>tol); 

%if relative tolerance instead...

ra = sum(S>max(S)*tol); 

%ra2 = sum(N>max(N)*tol);

%invert matrix containing singular values, up to a tolerance



if visualise
   SS = diag(S);
   % Sv = vecnorm(SS*V',2,1);

    figure(55)
   % clf;
    semilogy(S,'o-');
    hold on
    semilogy(ra*ones(1,2),logspace(-15,5,2),'r--')
    c = max(S)/min(S);  %Condition number
    str = sprintf('Self condition number %1.3e, max sing %1.3e',c,max(S));
    title(str,'interpreter','latex');
    grid on
    xlabel('$j$','interpreter','latex')
    ylabel('$\sigma_j$','interpreter','latex')

    % figure(56)
    % SV = SS*V';
    % semilogy(abs(SV(:,end-50:end)'))
    % 
    % 
    % figure(56)
    % semilogy(abs(UU(end-3:end,:)'))

    % figure(57)
    % semilogy(abs(diff(S)));
    % hold on
    % title('Decay rate of sing vals','interpreter','latex')
    % 
    % 
    % figure(58)
    % semilogy(abs(diff(S)./S(2:end)));
    % hold on
    % title('Relative decay rate of sing vals','interpreter','latex')
end

if truncate_last
    ra = ra-1; 
end
S = S(1:ra);  
iS = 1./S; % rank
Y = V(:,1:ra)*diag(iS); 
Bi1 = UU(:,1:ra); 

%Check if obeys with convergence predicution for the singular values
% figure()
% p = polyfit(sqrt(1:ra-1),log(S(1:end-1)),1);
% c2 = p(1)/log(0.59)
% y = exp(p(2))*exp(p(1)*sqrt(1:ra-1));
% 
% figure()
% % plot(log(S(1:end-1)));
% % hold on
% % plot(log(y))
% semilogy(S(1:end-1))
% hold on
% semilogy(y)
%y = exp(log(0.59)*2*sqrt(1:ra-1))

%semilogy(y)

end