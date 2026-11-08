function system = bubble_explicit_system(geometry,opt,Fg)
%BUBBLE_EXPLICIT_SYSTEM Cartesian cap rows and independent [lambda;U].
% Rows: cap velocity(3Mc), mobile normal(Mm), mobile traction(2Mm), force(3).
% FH=-KN'*lambda, so the last right-hand-side block is +Fg.
if nargin==0
    self_test();
    system = [];
    return
end
Fg = Fg(:);
assert(numel(Fg)==3 && isreal(Fg) && all(isfinite(Fg)) && ...
    norm(Fg)>0 && norm(Fg(1:2))<1e-13*norm(Fg), ...
    'Prescribed force must be nonzero and axial.');
assert(isequal(opt.refinementFamilies,"S") && ...
    opt.derivative_direction=="none" && ...
    geometry.coefficientCount==3*geometry.sourceCount, ...
    'bubble_explicit_system:StokesletsOnly', ...
    'The explicit bubble formulation requires Stokeslets only.');
p = bubble_boundary_operators(geometry,opt);
M = geometry.collocationCount;
N = geometry.sourceCount;
Mc = geometry.capCount;
Mm = geometry.mobileCount;
rigidMap = getKmat(geometry.sourcePoints,opt.centre);
KN = rigidMap(:,1:3);
Bm = p.Nm*p.KMm;
system.matrix = [p.Gc,-p.KMc; ...
    p.Nm*p.G(p.mobileVectorRows,:),-Bm; ...
    p.Dm*p.Tm,zeros(2*Mm,3); ...
    KN',zeros(3)];
system.rhs = [zeros(3*M,1);Fg];
system.forceMap = KN';
system.momentMap = rigidMap';
system.basis = p.basis;
system.capVelocityRows = 1:3*Mc;
system.mobileNormalRows = 3*Mc+(1:Mm);
system.mobileTractionRows = 3*Mc+Mm+(1:2*Mm);
system.forceRows = 3*M+(1:3);
system.blockCounts = [3*Mc,Mm,2*Mm,3];
system.Bm = Bm;
system.sourceCount = N;
system.coefficientCount = 3*N;
system.collocationCount = M;
assert(isequal(size(system.matrix),[3*M+3,3*N+3]));
assert(isequal(size(system.rhs),[3*M+3,1]) && sum(system.blockCounts)==3*M+3);
assert(norm(Bm-geometry.normals(geometry.mobile,:),'fro')<1e-12);
end

function self_test()
% Check cap equivalence, dimensions, force sign and endpoint empty blocks.
for levels = [0 1]
    for phi = [0 pi/3 pi]
        if levels>0 && (phi==0 || phi==pi), continue; end
        o = bubble_options(struct('sourceRequest',12,'phi',phi, ...
            'n_level',levels,'n_az_source',6,'n_az_colloc',8));
        g = bubble_geometry(o);
        s = bubble_explicit_system(g,o,[0;0;-1]);
        old = bubble_system(g,o,[0;0;-1]);
        alpha = sin((1:g.coefficientCount).');
        lambda = alpha-old.forceCompletion*(old.forceMap*alpha)+old.lambda0;
        U = -o.completionScale*old.forceMap*alpha;
        residual = s.matrix*[lambda;U]-s.rhs;
        previous = old.matrix*alpha-old.rhs;
        assert(abs(norm(residual(1:end-3))-norm(previous))< ...
            1e-11*max(1,norm(previous)));
        capResidual = reshape(residual(s.capVelocityRows),3,[]).';
        p = bubble_boundary_operators(g,o);
        assert(norm(sum(capResidual.*g.normals(g.cap,:),2)- ...
            previous(find(g.cap)))<1e-11); %#ok<FNDSB>
        capTangents = [sum(capResidual.*p.tangentOne(g.cap,:),2), ...
            sum(capResidual.*p.tangentTwo(g.cap,:),2)];
        assert(norm(reshape(capTangents.',[],1)- ...
            previous(g.collocationCount+(1:2*g.capCount)))<1e-11);
        assert(norm(residual(s.forceRows))<1e-12);
        expected = sum(reshape(lambda,3,[]),2);
        assert(norm(s.forceMap*lambda-expected)<1e-12);
        assert(all(s.matrix([s.mobileTractionRows,s.forceRows],end-2:end)==0,'all'));
    end
end
% An enclosing spherical design integrates a centred Stokeslet's force exactly.
[unit,~] = stagnant_cap_design(50);
f = [0.3;-0.4;1];
t = reshape(drop_traction_block([0 0 0],unit,unit)*f,3,[]).';
FH = (4*pi/size(unit,1))*(t.'*ones(size(unit,1),1));
assert(norm(FH+f)<1e-11);
Q = [0 0 1;1 0 0;0 1 0];
r = [1;2;3];
assert(abs(norm(Q*r)^2-norm(r)^2)<1e-13);
assert(abs(norm(diag([2 1 1])*Q*r)^2-norm(r)^2)>1);
fprintf('bubble_explicit_system self-test passed.\n');
end
