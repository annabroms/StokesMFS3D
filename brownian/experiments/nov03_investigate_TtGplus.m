close all; clear;
% This script aims at investigating the matrix T^TG^+. In principle, this
% approximates the inverse of the single layer S, mapping outer surface to
% outer surface. However, T^TG^+ is rank-deficient due to the
% rectangularity of G and therefore not a reliable approximation

%% Set initial particle coordinates
P=1;
delta = 5; 
if P==1
    q = [0 0 0];
else
    delta = 3; %initial particle particle distance
    q = [0 0 0; 2+delta 0 0];
end

%% Discretize!
%Rp = 0.33;
Rp = 0.15; %Proxy sphere radius -- very coarse resolution
N = 50;  % Number of proxy sources

%Rp = 0.10;
% N = 20; 

% Rp = 0.63; %fine resolution
% N = 700;
a = 1.5; %Determines oversampling factor for the collocation points
[rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a); %Assign source and collocation points
N = size(rvec_in,1)/P;
M = size(rvec_out,1)/P; 
%rvec_in = rvec_in(1,:);
n = repmat(rvec_out(1:M,:),P,1); %set surface normals

%% Get layer potential  and traction and their pseudoinverses
Tt = getTractionMat(rvec_in,rvec_out,n); % In the notes, we denote this matrix T'.
%T = getStresslet(rvec_in,rvec_out, repmat(rvec_in(1:N,:),P,1));
%This transpose is the same thing!

G = stokes_SLP_mat(rvec_in, rvec_out);
visualise = 1; 
tol = 1e-15;

[Y,U]  = getPseudoFactors(G,tol,visualise);
[YT,UT]  = getPseudoFactors(-Tt',tol,visualise);

[UU,SS,V] = svd(G);
Tprime = UU(1:end-1)*UU(1:end-1)'*Tt*V(:,1:end-1)*V(:,1:end-1)';
% atempt at projecting T' on the singular vectors of G

%% Check pseudoinverse of G
reg_err1 = norm(Y*(U'*G)-eye(3*N*P)); % small
reg_err2 = norm(G*Y*U'-eye(3*M*P)); % not small, not approx of identity
figure()
subplot(1,2,1)
imagesc(log10(abs(Y*(U'*G))))
colorbar
title('G^+G Should be good approximation to identity')
subplot(1,2,2)
imagesc(log10(abs(G*Y*U')))
colorbar
title('GG^+')

%% Solve resistance problem using G^+ only, or by first computing tractions. Same result. 
Kout = getKmat(rvec_out,q);
Kin = getKmat(rvec_in,q);
Urand = rand(6*P,1); %random translational and angular velocities
F_orig = Kin'*Y*(U'*Kout*Urand);
F_traction = 4*pi*Kout'*Tt*Y*(U'*Kout*Urand)/M;
%% Compare matrices

S_RPY = generate_RPY_matrix(rvec_out,0.03);
%B2 = (4*pi/M)*S*YT*UT';
%B2 = -S*YT*UT';
S = stokes_SLP_mat(rvec_out, rvec_out);

% Set self interaction terms in S to zero (they blow up)
for i = 1:M
    S(3*i-2:3*i, 3*i-2:3*i) = 0; %Binv(3*i-2:3*i, 3*i-2:3*i);
end

%Should use this scaling here too?
B = -(4*pi/M)*(Tt*Y)*U'; % T'G^+
[YB,UB]  = getPseudoFactors(B,tol,1);
% [V,D] = eig(B);
% e = diag(D); 

%build B step by step (as order of multiplication matters for pseudoinv)
% ek = zeros(3*M,1);
% for k = 1:3*M
%     ek(:) = 0; 
%     ek(k) = 1; 
%     B2(:,k) = -(4*pi/M)*T'*Y*(U'*ek);
% end

% Same pseudoinverse test as above
% figure()
% imagesc(log10(abs(YB*(UB'*B))))
% imagesc(log10(abs(B*YB*UB')))

% [UB,SB,VB] = svd(B); 
% r = 2058;
% Breg = UB(:,1:r)*SB(1:r,1:r)*VB(:,1:r)';

Binv = YB*UB';
%Binv = B\eye(M*P*3); % rank deficient...
%% Visualise inv of T^TG+ and compare against S
figure()
imagesc(log10(abs(Binv)))
colorbar
title('inv of T^TG+')
clim([-9,0])

figure()
imagesc(log10(abs(S))) %looks the same, but scaling is a bit off?
colorbar
title('S')
clim([-9,0])

%diff_mat = (S-Binv)./S; %off by some scaling?
diff_mat = (S-Binv);

figure()
imagesc(log10(abs(diff_mat)))
colorbar
title('Difference between S and pseudoinverse of T^TG^+: not matching and off by scaling')

%% Attempt at building RBM correlation expression  
K  = getKmat(rvec_in,q);
KGplus = K'*Y*U';
M_corr = M*KGplus*Binv*KGplus'*M'; 







