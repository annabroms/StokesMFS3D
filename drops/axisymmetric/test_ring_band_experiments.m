function test_ring_band_experiments()
%TEST_RING_BAND_EXPERIMENTS Check band mapping and small shared-grid solves.
% No inputs or outputs. Uses one sphere, unequal ring populations, Cartesian
% node-interleaved coefficients, and normalized independent boundary checks.
opt = struct('a',1,'centre',[0 0 0],'phi',pi/4,'mu',1,'Fg',[0;0;-1], ...
    'proxyRadius',0.45,'Nu',6,'Nv',4,'sourceWidth',0.4,'sourceSpacing',0.25, ...
    'collocWidth',0.5,'collocSpacing',0.3,'Mreq',60,'collocationLevels',2, ...
    'collocationSigma',2,'collocationAzimuths',12,'ringMode',"full", ...
    'validationMultiplier',5,'validationOffset',0.01,'auditFactor',2);
depths = [0.15 0.2];
d = ring_band_setup(opt,depths,{"SS_y"});
b = d.band;
assert(numel(unique(b.sourceCounts))>1);
[basis,C,w,~,H] = ring_mixed_basis(d.Y0,d.id0,b.sources,b.sourceRing, ...
    opt.centre,opt.a,["S","S_y"],"full");
R0 = max(d.id0); R = numel(b.sourceCounts);
coeff = zeros(size(C,2),1);
coeff(3:3:3*R0) = (1:R0)'/R0;
coeff(3*R0+3:3:3*R0+3*R) = (1:R)'/R;
expectedForce = sum(accumarray(d.id0,1).*(1:R0)'/R0) + ...
    sum(b.sourceCounts.*(1:R)'/R);
assert(abs(w.'*coeff-expectedForce)<1e-12);
coeff(:) = 0;
coeff(3*R0+3*R+(2:3:3*R)) = 1;
lambda = reshape(C*coeff,3,[]);
N0 = size(d.Y0,1); N = size(b.sources,1);
q = basis(3).normalDirections;
expectedTorque = sum(cross(q,lambda(:,N0+N+(1:N)).',2),1).'*b.depth;
assert(norm(H(1:3,:)*coeff)<1e-12);
assert(norm(H(4:6,:)*coeff-expectedTorque)<1e-12 && norm(expectedTorque)>1e-3);

other = ring_band_geometry(opt.centre,opt.a,opt.phi,depths(2), ...
    opt.sourceWidth,opt.sourceSpacing,opt.collocWidth,opt.collocSpacing);
scaled = opt.centre+(opt.a-depths(2))/(opt.a-b.depth)*(b.sources-opt.centre);
assert(norm(other.sources-scaled,'fro')<1e-12 && isequal(other.sourceRing,b.sourceRing));
for mode = ["full","swirl_free"]
    r = ring_family_case(d.Y0,d.id0,other.sources,b.sourceRing,opt.centre,opt.a, ...
        opt.mu,opt.Fg,"SS_y",d.cg,d.op,d.val,d.audit,d.Uexact, ...
        1/(6*pi),1e-10,mode,opt.phi);
    assert(all(isfinite([r.U;r.residual(:);r.auditResidual(:)])) && r.forceResidual<1e-12);
    for k = 1:3
        assert(abs(max(r.mapResiduals{k})-r.auditResidual(k))<1e-13);
    end
    idx = d.audit.cap;
    theta = d.audit.theta(idx);
    near = abs(theta-opt.phi)<=0.1;
    assert(any(near));
    assert(abs(max(r.mapResiduals{1}(near))-r.regional(1,1))<1e-13);
    if mode=="swirl_free"
        r2 = ring_family_case(d.Y0,d.id0,scaled,b.sourceRing,opt.centre,opt.a, ...
            opt.mu,opt.Fg,"SS_y",d.cg,d.op,d.val,d.audit,d.Uexact, ...
            1/(6*pi),1e-10,mode,opt.phi);
        assert(norm(r.U-r2.U)<1e-12 && norm(r.auditResidual-r2.auditResidual)<1e-12);
    end
end
fprintf('test_ring_band_experiments passed.\n');
end
