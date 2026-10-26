function A = stokes_directional_stokeslet_matrix(y,x,q,d,kind,n)
%STOKES_DIRECTIONAL_STOKESLET_MATRIX Source derivatives of Stokeslet fields.
% The viscosity-independent, force-normalized blocks include the per-source
% distance d after differentiation. Rows/columns are node-interleaved.
if nargin==0
    self_test(); A = []; return
end
if nargin<6, n = []; end
N = size(y,1); M = size(x,1);
assert(isequal(size(y),[N 3]) && isequal(size(x),[M 3]) && ...
    isequal(size(q),[N 3]) && isequal(size(d),[N 1]) && ...
    all(isfinite(y),'all') && all(isfinite(x),'all') && ...
    all(isfinite(q),'all') && all(isfinite(d)) && all(d>0));
assert(max(abs(vecnorm(q,2,2)-1))<1e-10);
kind = string(kind);
if kind=="traction"
    assert(isequal(size(n),[M 3]) && all(isfinite(n),'all'));
elseif ~ismember(kind,["velocity","pressure"])
    error('stokes_directional_stokeslet_matrix:BadKind', ...
        'Use velocity, pressure, or traction.');
end
r = cell(1,3);
for k = 1:3
    r{k} = x(:,k)-y(:,k).';
end
r2 = r{1}.^2+r{2}.^2+r{3}.^2;
if any(r2==0,'all')
    error('stokes_directional_stokeslet_matrix:CoincidentPoints', ...
        'A source coincides with a target.');
end
ir3 = r2.^(-3/2);
ir5 = ir3./r2;
rq = r{1}.*q(:,1).'+r{2}.*q(:,2).'+r{3}.*q(:,3).';
if kind=="pressure"
    A = zeros(M,3*N);
    for k = 1:3
        A(:,k:3:end) = d.'.*(-q(:,k).'.*ir3+3*r{k}.*rq.*ir5)/(4*pi);
    end
    return
end
A = zeros(3*M,3*N);
if kind=="traction"
    rn = r{1}.*n(:,1)+r{2}.*n(:,2)+r{3}.*n(:,3);
    qn = n*q.';
    ir7 = ir5./r2;
end
for i = 1:3
    for k = 1:3
        if kind=="velocity"
            value = ((i==k)*rq-q(:,i)'.*r{k}-r{i}.*q(:,k)').*ir3 ...
                +3*r{i}.*r{k}.*rq.*ir5;
            value = value/(8*pi);
        else
            value = 3/(4*pi)*(qn.*r{i}.*r{k}+ ...
                rn.*(q(:,i)'.*r{k}+r{i}.*q(:,k)')).*ir5 ...
                -15/(4*pi)*rn.*r{i}.*r{k}.*rq.*ir7;
        end
        A(i:3:end,k:3:end) = d.'.*value;
    end
end
end

function self_test()
%SELF_TEST Check source signs, pressure, traction, and finite-difference order.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'common','dense'),fullfile(root,'common','eval'), ...
    fullfile(root,'drops'));
x = [1.2 0.4 0.7;0.3 -0.8 1.1;-0.7 0.6 -1.2];
n = x./vecnorm(x,2,2);
b = [0.4;-0.8;0.3];
directions = [[1 0 0];[0.3 0.4 sqrt(0.75)]];
sources = [[0.6 0 0];[1.35 0 0]];
d = 0.2;
steps = 10.^(-2:-1:-6);
for side = 1:2
    y = sources(side,:);
    q = directions(side,:);
    exactU = stokes_directional_stokeslet_matrix(y,x,q,d,"velocity")*b;
    exactP = stokes_directional_stokeslet_matrix(y,x,q,d,"pressure")*b;
    exactT = stokes_directional_stokeslet_matrix(y,x,q,d,"traction",n)*b;
    errorU = zeros(size(steps)); errorP = errorU; errorT = errorU;
    for j = 1:numel(steps)
        h = steps(j);
        yp = y+h*q; ym = y-h*q;
        assert((norm(yp)<1 && norm(ym)<1) || ...
            (norm(yp)>1 && norm(ym)>1));
        fdU = d*(stokes_SLP_mat(yp,x)-stokes_SLP_mat(ym,x))*b/(2*h);
        fdP = d*(getPressure(b,yp,x)-getPressure(b,ym,x))/(2*h);
        fdT = d*(drop_traction_block(yp,x,n)- ...
            drop_traction_block(ym,x,n))*b/(2*h);
        errorU(j) = norm(fdU-exactU,inf);
        errorP(j) = norm(fdP-exactP,inf);
        errorT(j) = norm(fdT-exactT,inf);
    end
    assert(errorU(3)<errorU(1)/100 && errorP(3)<errorP(1)/100 && ...
        errorT(3)<errorT(1)/100);
    assert(max([errorU(4),errorP(4),errorT(4)])<2e-8);
end
fprintf('stokes_directional_stokeslet_matrix self-test passed.\n');
end
