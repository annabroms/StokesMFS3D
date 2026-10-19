function [u,p,t] = drop_image_field(x,n,region,grid,lambda,gamma,opt)
%DROP_IMAGE_FIELD Evaluate exterior or one drop's interior flow representation.
% region=0 selects global exterior. region=i uses outer sources of drop i
% plus incident exterior fields from all other drops. Pressure tends to zero
% at infinity in every free-space representation.
if nargin == 0
    self_test(); u=[]; p=[]; t=[]; return
end
P = numel(grid.body);
if isempty(n), n = zeros(size(x)); end
if ~isequal(size(x),size(n)) || region<0 || region>P
    error('drop_image_field:BadInput','Bad target normals or region.');
end
u = zeros(size(x)); p = zeros(size(x,1),1); t = zeros(size(x));
for j = 1:P
    if j==region
        basis = grid.body(j).interiorBasis;
        z = gamma{j};
        velocityFactor = 1/opt.interiorViscosity(j);
        pressureFactor = 1;
        tractionFactor = 1;
    else
        basis = grid.body(j).basis;
        z = lambda{j};
        velocityFactor = 1;
        pressureFactor = 1;
        tractionFactor = 1;
        if region>0
            pressureFactor = opt.interiorViscosity(region);
            tractionFactor = pressureFactor;
        end
    end
    coeff = image_source_coefficients(basis,z);
    one.body.basis = basis;
    for first = 1:opt.batchSize:size(x,1)
        rows = first:min(first+opt.batchSize-1,size(x,1));
        xx = x(rows,:);
        u(rows,:) = u(rows,:)+velocityFactor* ...
            direct_image_velocity(basis,coeff,xx);
        p(rows) = p(rows)+pressureFactor*image_pressure(one,{z},xx,1);
        if nargout>2
            t(rows,:) = t(rows,:)+tractionFactor* ...
                image_traction(one,{z},xx,n(rows,:),1);
        end
    end
end
end

function self_test()
% Verify the base Stokeslet pressure and viscosity-scaled incident field.
q = [0 0 0;3 0 0]; a = ones(2,1);
opt = drop_image_options(struct('N',30,'useImages',false,'batchSize',2),2,a);
g = drop_image_discretization(q,a,opt);
lam = cell(2,1); gam = cell(2,1);
for i = 1:2
    lam{i} = zeros(3*size(g.body(i).basis(1).points,1),1);
    gam{i} = zeros(3*size(g.body(i).interiorBasis(1).points,1),1);
end
lam{2}(1)=1;
x = [0 0 1;0 1 0]; n = x;
[u0,p0] = drop_image_field(x,n,0,g,lam,gam,opt);
[ui,interiorPressure] = drop_image_field(x,n,1,g,lam,gam,opt);
assert(norm(u0-ui,inf)<1e-13);
assert(norm(interiorPressure-opt.interiorViscosity(1)*p0,inf)<1e-13);
fprintf('drop_image_field self-test passed.\n');
end
