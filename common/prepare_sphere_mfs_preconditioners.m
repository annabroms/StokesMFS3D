function cache = prepare_sphere_mfs_preconditioners(rin0,rout0)
%PREPARE_SPHERE_MFS_PRECONDITIONERS Factor the image-free one-body maps once.
%   cache = prepare_sphere_mfs_preconditioners(rin0,rout0)
% rin0 and rout0 are N-by-3 and M-by-3 body-frame nodes about the origin.
% cache.resistance and cache.mobility may be passed as opt.precond to the
% corresponding basic MFS solvers for any translation or rigid rotation of
% these same nodes. The truncation threshold matches oneBodyPrecondRes/Mob.

validateattributes(rin0,{'numeric'},{'real','finite','2d','ncols',3});
validateattributes(rout0,{'numeric'},{'real','finite','2d','ncols',3});
N = size(rin0,1);
M = size(rout0,1);
assert(N >= 3 && M >= N,'Expected at least three sources and M >= N.');

setupTimer = tic;
Sblock = stokes_SLP_mat(rin0,rout0);
Kin = getKmat(rin0,[0,0,0]);
Kout = getKmat(rout0,[0,0,0]);
[Yres,Ures,singularValues] = getPseudoFactors(Sblock,1e-15,0);
Ures = Ures.';

LL = Kin*((Kin.'*Kin)\Kin.');
Amob = Sblock*(eye(3*N)-LL)+Kout*Kin.';
[Ymob,Umob] = getPseudoFactors(Amob,1e-15,0);
Umob = Umob.';

common = struct('rin0',rin0,'rout0',rout0,'Sblock',Sblock);
cache.resistance = common;
cache.resistance.mode = 'resistance';
cache.resistance.Y = Yres;
cache.resistance.UU = Ures;
cache.resistance.singularValues = singularValues;
cache.mobility = common;
cache.mobility.mode = 'mobility';
cache.mobility.Y = Ymob;
cache.mobility.UU = Umob;
cache.mobility.LL = LL;
cache.mobility.Kin = Kin;
cache.setupSeconds = toc(setupTimer);
cache.sourceCount = N;
cache.collocationCount = M;
cache.svdCutoff = 1e-15;
end
