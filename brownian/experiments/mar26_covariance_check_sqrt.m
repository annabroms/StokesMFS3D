% Compare empirical covariance of Lanczos samples against the dense C matrix.
% The Krylov solve uses the matrix-free dense = 0 path, while C is assembled
% once for post-processing and covariance comparisons.

clear;
close all;


rng(5);

sphere = 0;
tdesign = 0;
dense = 0;
fast = 1;
fmm = 0;

nruns = 100; % Increase for a smoother covariance estimate.
tol = 1e-8;
maxit = 500;
trunc = 1e-12;

P = 2;
delta = 5;
E0 = [0.4 0.6 1.0];
sep = 0.125;

if sphere
    E0 = [1 1 1];
    sep = 0.3;
end

Nv = 15;
a = 2;

rng(8);
[~, R, qvec, ~] = ellipsoid_cluster(E0,P,delta);
disp('Geometry created')

for k = 1:P
    R{k} = eye(3); 
end

N1 = 0.75*Nv;
N2 = Nv;
if sphere && tdesign
    Nt = 50;
    Rp = 0.15;
    q = cell2mat(cellfun(@transpose, qvec, 'UniformOutput', false).');
    [rvec_in,rvec_out,~] = init_spheres(q,Rp,Nt,a);
    M = size(rvec_out,1)/P;
    n_out = rvec_out - kron(q,ones(M,1));
    w = (4*pi/M)*ones(P*M,1);
    Rprecond = [];
else
    if sphere
        N1 = 6;
        N2 = 8;
        sep = 0.85;
    end
    [rvec_in,rvec_out,~,~,~,w,n_out] = getEllipsoidGrids(E0,P,delta,N1,N2,sep,R,qvec);
    Rprecond = R;
end

M = size(rvec_out,1)/P;
N = size(rvec_in,1)/P;
ndof = 3*N*P;

figure(1)
scatter3(rvec_in(:,1),rvec_in(:,2),rvec_in(:,3),'filled')
hold on
scatter3(rvec_out(:,1),rvec_out(:,2),rvec_out(:,3),'filled')
axis equal
title('Input and collocation nodes')

% Dense reference matrix used only for covariance checks.
[C,~,~] = getSymTWS(rvec_out,rvec_in,n_out,w);

% One-body preconditioner.
Ci = getSymTWS(rvec_out(1:M,:),rvec_in(1:N,:),n_out(1:M,:),w(1:M));
[Vi,Di] = eig(Ci);
di = diag(Di);

ind = find(real(di) > trunc);
diff_set = setdiff(1:numel(di),ind);
dsqrt = sqrt(di);
dsqrt_inv = 1./dsqrt;
dsqrt_inv(diff_set) = 0;
Pi = diag(dsqrt_inv)*Vi';

dsqrt_plus = dsqrt;
dsqrt_plus(diff_set) = 0;
Pi_plus = Vi*diag(dsqrt_plus);
PPplus = @(x) apply_block_diag(x, Pi_plus, P, Rprecond);

vars.fmm = fmm;
vars.eps = 1e-10;
Cfun = @(x) getPrecondTG(x,P,rvec_in,rvec_out,n_out,w,Pi,vars,Rprecond);

dW_hist = randn(ndof,nruns);
y_hist = zeros(ndof,nruns);
iter_err_hist = cell(nruns,1);
krylov_iters = zeros(nruns,1);
cov_rel_err = zeros(nruns,1);

norm_type = 2;
for ir = 1:nruns
    dW = dW_hist(:,ir);

    if fast
        [y,iter_err,iters] = lanczosSqrt_fast(Cfun,dW,tol,maxit,PPplus);
    else
        [y,iter_err,iters] = lanczosSqrt(Cfun,dW,tol,maxit,PPplus);
    end

    y_hist(:,ir) = y;
    iter_err_hist{ir} = iter_err(:);
    krylov_iters(ir) = iters;

    Yrun = y_hist(:,1:ir);
    C_emp_run = (Yrun*Yrun.')/ir; % transpose, not conjugate transpose
    cov_rel_err(ir) = norm(C_emp_run - C,norm_type)/max(norm(C,norm_type),eps);
end

y_mean = mean(y_hist,2);
C_emp = (y_hist*y_hist.')/nruns;
Y_centered = y_hist - y_mean;
C_emp_centered = (Y_centered*Y_centered.')/max(nruns-1,1);


rel_cov_err = norm(C_emp - C,norm_type)/max(norm(C,norm_type),eps);
rel_cov_err_centered = norm(C_emp_centered - C,norm_type)/max(norm(C,norm_type),eps);
rel_diag_err = norm(diag(C_emp) - diag(C))/max(norm(diag(C)),eps);
mean_rel = norm(y_mean)/max(mean(vecnorm(y_hist)),eps);

fprintf('mar26_covariance_check_sqrt:\n');
fprintf('  ndof = %d\n', ndof);
fprintf('  nruns = %d\n', nruns);
fprintf('  mean Krylov iterations = %.2f\n', mean(krylov_iters));
fprintf('  max Krylov iterations  = %d\n', max(krylov_iters));
fprintf('  rel error, raw second moment      = %.3e\n', rel_cov_err);
fprintf('  rel error, centered second moment = %.3e\n', rel_cov_err_centered);
fprintf('  rel diagonal error                          = %.3e\n', rel_diag_err);
fprintf('  relative sample mean norm                   = %.3e\n', mean_rel);
if ~isreal(y_hist)
    fprintf('  Samples are complex; moments use transpose rather than conjugate transpose.\n');
end

figure(2)
subplot(2,2,1)
imagesc(real(C))
axis image
colorbar
title('real(C)')

subplot(2,2,2)
imagesc(real(C_emp))
axis image
colorbar
title(sprintf('real(C_{emp}), nruns = %d', nruns))

subplot(2,2,3)
imagesc(log10(abs(C_emp - C) + eps))
axis image
colorbar
title('log10 |C_{emp} - C|')

subplot(2,2,4)
semilogy(cov_rel_err,'LineWidth',1.5)
grid on
xlabel('Sample index')
ylabel('Rel. error')
title('Running covariance error')

figure(3)
plot(real(diag(C)),real(diag(C_emp)),'.','MarkerSize',12)
hold on
plot(xlim,xlim,'k--')
grid on
axis equal
xlabel('diag(C)')
ylabel('diag(C_{emp})')
title('Diagonal comparison')

results.C = C;
results.C_emp = C_emp;
results.C_emp_centered = C_emp_centered;
results.y_hist = y_hist;
results.dW_hist = dW_hist;
results.iter_err_hist = iter_err_hist;
results.krylov_iters = krylov_iters;
results.cov_rel_err = cov_rel_err;
results.rel_cov_err = rel_cov_err;
results.rel_cov_err_centered = rel_cov_err_centered;
results.rel_diag_err = rel_diag_err;
results.y_mean = y_mean;
