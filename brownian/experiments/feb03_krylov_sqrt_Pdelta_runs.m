% Krylov sqrt iterations vs P and delta (fixed tolerance), one-body precond
% Lanczos

clear; close all;

rng(5);

Pvec = 5:5:40;
deltavec = [0.5 1 5 10];
Nruns = 5;

fmm = 0; % only activate if many particles

% Fixed resolution and parameters (same baseline as jan28_test_krylov_sqrt_Plarge)
Rp = 0.30;
N = 100;
Rp = 0.15;
N = 50; 
a = 2;

tol = 1e-8;
tol = 1e-12;
maxit = 500;
err_target = 1e-4;

vars.fmm = fmm;
vars.eps = 1e-10; % for use in FMM

nP = numel(Pvec);
nD = numel(deltavec);
iter_req = nan(nP, nD, Nruns);

for ip = 1:nP
    P = Pvec(ip);
    for id = 1:nD
        delta = deltavec(id);

        for ir = 1:Nruns
            ir
            if P == 2
                q = [0 0 0; 2+delta 0 0];
            elseif P == 1
                q = [0 0 0];
            else
                [q,~] = grow_cluster(P,delta);
            end

            [rvec_in,rvec_out,opt] = init_spheres(q,Rp,N,a);
            M = size(rvec_out,1)/P;
            N = size(rvec_in,1)/P;
            n_out = repmat(rvec_out(1:M,:),P,1);

            dW2 = rand(3*N*P,1);

            % Block-diagonal preconditioner
            Ti = stokes_DLP(rvec_out(1:M,:),rvec_in(1:N,:),rvec_out(1:M,:));
            Si = generate_stokes_mat(rvec_in(1:N,:), rvec_out(1:M,:));
            Ai = Ti*Si;
            Bi = (Ai+Ai')/2;
            [Ve,De] = eig(Bi);
            d = diag(De);

            trunc = tol; % fixed truncation at same level as Krylov tolerance
            ind = find(real(d)>trunc);
            diff_set = setdiff(1:length(d),ind);
            dsqrt = 1./sqrt(d);
            dsqrt(diff_set) = 0;
            Ci = diag(dsqrt)*Ve';

            dsqrt_plus = sqrt(d);
            dsqrt_plus(diff_set) = 0;
            Ci_plus = Ve*diag(dsqrt_plus);
            Cplus = @(x) apply_block_diag(x, Ci_plus, P);

            CBC = @(x) getPrecondTG(x,P,rvec_in,rvec_out,n_out,Ci,vars);

            [~,iter_err] = KrylovSqrtMsing(CBC,dW2,tol,maxit,Cplus);

            idx = find(iter_err <= err_target, 1, 'first');
            if isempty(idx)
                iter_req(ip,id,ir) = NaN; % did not reach target tolerance
            else
                iter_req(ip,id,ir) = idx;
            end
        end
    end
end

function y = apply_block_diag(x, B, P)
%APPLY_BLOCK_DIAG Apply block-diagonal operator with P identical blocks B.
% x is (3N*P) x 1, B is (3N) x (3N).
N = size(B,1);
y = zeros(size(x));
for k = 1:P
    idx = (k-1)*N + (1:N);
    y(idx) = B * x(idx);
end
end
% Visualize: median iterations to reach target + uncertainty (IQR)
med_iters = mean(iter_req, 3);
q25 = prctile(iter_req, 25, 3);
q75 = prctile(iter_req, 75, 3);
iqr_iters = q75 - q25;

%%

figure(1); clf

% Sizes
[nP, nD] = size(med_iters);

% Helper: build bin edges from centers (works for nonuniform spacing too)
centers_to_edges = @(c) [ ...
    c(1) - 0.5*(c(2)-c(1)), ...
    0.5*(c(1:end-1) + c(2:end)), ...
    c(end) + 0.5*(c(end)-c(end-1)) ];

% Decide whether deltavec/Pvec are edges or centers
if numel(deltavec) == nD+1 && numel(Pvec) == nP+1
    % --- Case A: already edges ---
    delta_edges = deltavec(:).';   % 1 x (nD+1)
    P_edges     = Pvec(:).';       % 1 x (nP+1)
elseif numel(deltavec) == nD && numel(Pvec) == nP
    % --- Case B: given centers; infer edges ---
    if nD < 2 || nP < 2
        error('Need at least 2 points in deltavec and Pvec to infer edges.');
    end
    delta_edges = centers_to_edges(deltavec(:).');
    P_edges     = centers_to_edges(Pvec(:).');
else
    error(['Dimension mismatch:\n' ...
           'size(med_iters) = [%d %d]\n' ...
           'numel(Pvec) = %d (expected %d or %d)\n' ...
           'numel(deltavec) = %d (expected %d or %d)\n'], ...
           nP, nD, numel(Pvec), nP, nP+1, numel(deltavec), nD, nD+1);
end

% Pad Z so surf has nP-by-nD patches
Zp = ones(nP+1, nD+1);
Zp(1:end-1, 1:end-1) = med_iters;

% Plot
surf(delta_edges, P_edges, Zp, 'EdgeColor','none');
view(0,90)
axis tight
xlabel('\delta')
ylabel('P')

cb = colorbar;
cb.Label.String = 'Mean iterations';
colormap(parula)

title(sprintf(['Krylov iterations to reach error %.1e (fixed eigval trunc = %.1e)\n' ...
    'Cells show median over %d runs; uncertainty is IQR = 75th–25th percentile'], ...
    err_target, tol, Nruns));

% Cell centers for text placement
dc = 0.5*(delta_edges(1:end-1) + delta_edges(2:end));  % 1 x nD
Pc = 0.5*(P_edges(1:end-1)     + P_edges(2:end));      % 1 x nP

% Text layer slightly above surface
zmax = max(med_iters(:), [], 'omitnan');
if isempty(zmax) || isnan(zmax), zmax = 1; end
zlevel = zmax * 1.01;

% Overlay values
for ip = 1:nP
    for id = 1:nD
        val    = med_iters(ip,id);
        spread = iqr_iters(ip,id);

        if ~isnan(val)
            txt = sprintf('%.0f\n(IQR %.0f)', val, spread);

            if val > 0.5*zmax
                tcolor = 'w';
            else
                tcolor = 'k';
            end

            text(dc(id), Pc(ip), zlevel, txt, ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'Color',tcolor, ...
                'FontWeight','bold', ...
                'Clipping','off');
        else
            text(dc(id), Pc(ip), zlevel, '—', ...
                'HorizontalAlignment','center', ...
                'VerticalAlignment','middle', ...
                'Color','k', ...
                'FontWeight','bold', ...
                'Clipping','off');
        end
    end
end
