function value = stokes_axial_quadrupole_field( ...
    sources,targets,axes,strengths,field,normals)
%STOKES_AXIAL_QUADRUPOLE_FIELD Sum trace-free axial Stokeslet quadrupoles.
% Q_s = (s.d_y)^2 G - Delta_y G/3 = (s.d_y)^2 G + D/3.
% The source direction s is unit length; strengths are physical N-by-3 rows.
if nargin == 0
    self_test(); value = []; return
end
if nargin < 6, normals = []; end
ns = size(sources,1);
assert(isequal(size(axes),[ns,3]) && ...
    isequal(size(strengths),[ns,3]));
assert(all(abs(vecnorm(axes,2,2)-1)<1e-10));
if field=="traction"
    assert(isequal(size(normals),size(targets)));
end

r = cell(1,3);
for j = 1:3
    r{j} = targets(:,j)-sources(:,j).';
end
r2 = r{1}.^2+r{2}.^2+r{3}.^2;
if any(r2==0,'all')
    error('stokes_axial_quadrupole_field:CoincidentPoints', ...
        'A target coincides with a source.');
end
ir3 = r2.^(-3/2);
ir5 = ir3./r2;
ir7 = ir5./r2;
a = r{1}.*axes(:,1).'+r{2}.*axes(:,2).'+r{3}.*axes(:,3).';
b = r{1}.*strengths(:,1).'+r{2}.*strengths(:,2).'+ ...
    r{3}.*strengths(:,3).';
h = sum(axes.*strengths,2).';
switch string(field)
    case "velocity"
        value = zeros(size(targets));
        for j = 1:3
            qss = strengths(:,j).'.*(-ir3+3*a.^2.*ir5) + ...
                2*axes(:,j).'.*h.*ir3 - ...
                6*a.*(axes(:,j).'.*b+r{j}.*h).*ir5 - ...
                3*r{j}.*b.*ir5+15*r{j}.*b.*a.^2.*ir7;
            dipole = (-strengths(:,j).'.*ir3+ ...
                3*r{j}.*b.*ir5)/(4*pi);
            value(:,j) = sum(qss/(8*pi)+dipole/3,2);
        end
    case "pressure"
        value = sum(((-6*h.*a-3*b).*ir5+ ...
            15*b.*a.^2.*ir7)/(4*pi),2);
    case "traction"
        value = zeros(size(targets));
        ir9 = ir7./r2;
        rn = r{1}.*normals(:,1)+r{2}.*normals(:,2)+ ...
            r{3}.*normals(:,3);
        sn = normals*axes.';
        cn = normals*strengths.';
        for j = 1:3
            qss = (2*axes(:,j).'.*h.*rn+ ...
                2*axes(:,j).'.*b.*sn+2*r{j}.*h.*sn).*ir5 - ...
                10*a.*(axes(:,j).'.*b.*rn+r{j}.*h.*rn+ ...
                r{j}.*b.*sn).*ir7 -5*r{j}.*b.*rn.*ir7+ ...
                35*r{j}.*b.*rn.*a.^2.*ir9;
            dipole = 3/(2*pi)*((strengths(:,j).'.*rn+ ...
                r{j}.*cn+b.*normals(:,j)).*ir5- ...
                5*b.*r{j}.*rn.*ir7);
            value(:,j) = sum(-3*qss/(4*pi)+dipole/3,2);
        end
    otherwise
        error('stokes_axial_quadrupole_field:BadField', ...
            'Expected velocity, pressure, or traction.');
end
end

function self_test()
% Check source derivatives, the Stokes PDE, and rotation covariance.
root = fileparts(fileparts(mfilename('fullpath')));
oldPath = path; cleanup = onCleanup(@() path(oldPath));
addpath(fullfile(root,'common','dense'),fullfile(root,'common','eval'));
y = [0.13 -0.21 0.08];
x = [1.1 0.5 -0.7;0.8 -0.9 0.6];
s = [1 2 3]/sqrt(14);
f = [0.4 -0.3 0.8];
n = x./vecnorm(x,2,2);
h = 1e-3;
b = struct('family',"S",'points',y, ...
    'normalDirections',zeros(0,3),'scale',1);
d = b; d.family = "D";
u = stokes_axial_quadrupole_field(y,x,s,f,"velocity");
p = stokes_axial_quadrupole_field(y,x,s,f,"pressure");
t = stokes_axial_quadrupole_field(y,x,s,f,"traction",n);
b.points = y+h*s;
up = reshape(stokes_SLP_mat(b.points,x)*f.',3,[]).';
pp = getPressure(f.',b.points,x);
tp = reshape(stokes_mfs_traction_matrix(b,x,n)*f.',3,[]).';
b.points = y-h*s;
um = reshape(stokes_SLP_mat(b.points,x)*f.',3,[]).';
pm = getPressure(f.',b.points,x);
tm = reshape(stokes_mfs_traction_matrix(b,x,n)*f.',3,[]).';
b.points = y;
u0 = reshape(stokes_SLP_mat(y,x)*f.',3,[]).';
p0 = getPressure(f.',y,x);
t0 = reshape(stokes_mfs_traction_matrix(b,x,n)*f.',3,[]).';
ud = reshape(stokes_mfs_velocity_matrix(d,x)*f.',3,[]).';
td = reshape(stokes_mfs_traction_matrix(d,x,n)*f.',3,[]).';
err = [norm(u-(up-2*u0+um)/h^2-ud/3,'fro'), ...
    norm(p-(pp-2*p0+pm)/h^2), ...
    norm(t-(tp-2*t0+tm)/h^2-td/3,'fro')];
% Rotate every geometric and strength vector by the same proper rotation.
angle = 0.47;
R = [cos(angle) -sin(angle) 0;sin(angle) cos(angle) 0;0 0 1];
uRot = stokes_axial_quadrupole_field(y*R.',x*R.',s*R.',f*R.',"velocity");
pRot = stokes_axial_quadrupole_field(y*R.',x*R.',s*R.',f*R.',"pressure");
tRot = stokes_axial_quadrupole_field( ...
    y*R.',x*R.',s*R.',f*R.',"traction",n*R.');
errRot = max([norm(uRot-u*R.','fro'),norm(pRot-p), ...
    norm(tRot-t*R.','fro')]);
if max(err)>1e-5 || errRot>1e-12
    error('stokes_axial_quadrupole_field:SelfTest', ...
        'Derivative or rotation check failed.');
end
fprintf('stokes_axial_quadrupole_field self-test: derivative %.3e, rotation %.3e.\n', ...
    max(err),errRot);
end
