function A = stokes_axial_quadrupole_matrix( ...
    sources,targets,axes,field,normals)
%STOKES_AXIAL_QUADRUPOLE_MATRIX Build dense trace-free axial quadrupole blocks.
% Velocity and traction rows are node-interleaved; pressure has one row per
% target. Columns are always node-interleaved Cartesian strengths.
if nargin == 0
    self_test(); A = []; return
end
if nargin < 5, normals = []; end
ns = size(sources,1);
nt = size(targets,1);
assert(isequal(size(axes),[ns,3]));
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
    error('stokes_axial_quadrupole_matrix:CoincidentPoints', ...
        'A target coincides with a source.');
end
ir3 = r2.^(-3/2);
ir5 = ir3./r2;
ir7 = ir5./r2;
a = r{1}.*axes(:,1).'+r{2}.*axes(:,2).'+r{3}.*axes(:,3).';
if field=="pressure"
    A = zeros(nt,3*ns);
    for k = 1:3
        b = r{k};
        h = axes(:,k).';
        A(:,k:3:3*ns) = ((-6*h.*a-3*b).*ir5+ ...
            15*b.*a.^2.*ir7)/(4*pi);
    end
    return
end
A = zeros(3*nt,3*ns);
if field=="traction"
    ir9 = ir7./r2;
    rn = r{1}.*normals(:,1)+r{2}.*normals(:,2)+ ...
        r{3}.*normals(:,3);
    sn = normals*axes.';
elseif field~="velocity"
    error('stokes_axial_quadrupole_matrix:BadField', ...
        'Expected velocity, pressure, or traction.');
end
for k = 1:3
    b = r{k};
    h = axes(:,k).';
    for j = 1:3
        if field=="velocity"
            qss = double(j==k)*(-ir3+3*a.^2.*ir5)+ ...
                2*axes(:,j).'.*h.*ir3 - ...
                6*a.*(axes(:,j).'.*b+r{j}.*h).*ir5 - ...
                3*r{j}.*b.*ir5+15*r{j}.*b.*a.^2.*ir7;
            dipole = (-double(j==k)*ir3+3*r{j}.*b.*ir5)/(4*pi);
            value = qss/(8*pi)+dipole/3;
        else
            qss = (2*axes(:,j).'.*h.*rn+ ...
                2*axes(:,j).'.*b.*sn+2*r{j}.*h.*sn).*ir5 - ...
                10*a.*(axes(:,j).'.*b.*rn+r{j}.*h.*rn+ ...
                r{j}.*b.*sn).*ir7 -5*r{j}.*b.*rn.*ir7+ ...
                35*r{j}.*b.*rn.*a.^2.*ir9;
            dipole = 3/(2*pi)*((double(j==k)*rn+ ...
                r{j}.*normals(:,k)+b.*normals(:,j)).*ir5- ...
                5*b.*r{j}.*rn.*ir7);
            value = -3*qss/(4*pi)+dipole/3;
        end
        A(j:3:3*nt,k:3:3*ns) = value;
    end
end
end

function self_test()
% Match dense velocity and traction actions to direct pair summation.
y = [0 0 0;0.2 -0.1 0.3];
x = [1.1 0.3 -0.2;-0.4 1.2 0.8];
s = [1 0 0;0 1 0];
n = x./vecnorm(x,2,2);
f = [0.4 -0.1 0.8;-0.3 0.7 0.2];
for field = ["velocity","pressure","traction"]
    A = stokes_axial_quadrupole_matrix(y,x,s,field,n);
    actual = A*reshape(f.',[],1);
    if field~="pressure", actual = reshape(actual,3,[]).'; end
    expected = stokes_axial_quadrupole_field(y,x,s,f,field,n);
    if norm(actual-expected,'fro')>1e-13
        error('stokes_axial_quadrupole_matrix:SelfTest', ...
            'Dense and direct quadrupole actions differ.');
    end
end
fprintf('stokes_axial_quadrupole_matrix self-test passed.\n');
end
